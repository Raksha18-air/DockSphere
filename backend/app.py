from flask import Flask, request, jsonify, send_from_directory, render_template, send_file
from flask_cors import CORS
from celery import Celery
from celery.result import AsyncResult
import os
import subprocess
import logging
import requests
import py3Dmol
import zipfile
import glob
from rdkit import Chem
from rdkit.Chem import AllChem
import uuid
from datetime import datetime
import shutil

app = Flask(__name__)
CORS(app)

# Configure logging
logging.basicConfig(level=logging.INFO)

# Celery configuration for async job processing
app.config['CELERY_BROKER_URL'] = os.environ.get('CELERY_BROKER_URL', 'redis://localhost:6379/0')
app.config['CELERY_RESULT_BACKEND'] = os.environ.get('CELERY_RESULT_BACKEND', 'redis://localhost:6379/0')

celery = Celery(app.name, broker=app.config['CELERY_BROKER_URL'])
celery.conf.update(app.config)

# Paths to required tools and directories
MGLTOOLS_PATH = os.environ.get("MGLTOOLS_PATH", "/opt/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24")
VINA_PATH = os.environ.get("VINA_PATH", "/usr/local/bin/vina")
MGL_PYTHON = os.environ.get("MGL_PYTHON", "/opt/mgltools/bin/python2.7")

BASE_FOLDER = os.environ.get("BASE_FOLDER", "/app/jobs")

# Create base directory
os.makedirs(BASE_FOLDER, exist_ok=True)

# ============================================
# CRITICAL FIX: Per-Job Isolated Directories
# ============================================
def create_job_directories(job_id):
    """
    Create isolated directories for each job to prevent file mixing
    """
    job_folder = os.path.join(BASE_FOLDER, job_id)
    
    upload_folder = os.path.join(job_folder, "uploads")
    protein_pdbqt_dir = os.path.join(job_folder, "proteins_pdbqt")
    ligand_pdbqt_dir = os.path.join(job_folder, "ligands_pdbqt")
    docking_results_dir = os.path.join(job_folder, "docking_results")
    
    for directory in [upload_folder, protein_pdbqt_dir, ligand_pdbqt_dir, docking_results_dir]:
        os.makedirs(directory, exist_ok=True)
    
    return {
        'job_folder': job_folder,
        'upload_folder': upload_folder,
        'protein_pdbqt_dir': protein_pdbqt_dir,
        'ligand_pdbqt_dir': ligand_pdbqt_dir,
        'docking_results_dir': docking_results_dir
    }

def cleanup_job_files(job_id, keep_results=True):
    """
    Clean up job files after completion
    keep_results: If True, only delete intermediate files
    """
    job_folder = os.path.join(BASE_FOLDER, job_id)
    
    if not os.path.exists(job_folder):
        return
    
    if keep_results:
        # Delete only uploads and intermediate PDBQT files
        for subdir in ['uploads', 'proteins_pdbqt', 'ligands_pdbqt']:
            path = os.path.join(job_folder, subdir)
            if os.path.exists(path):
                shutil.rmtree(path)
    else:
        # Delete entire job folder
        shutil.rmtree(job_folder)

#Function to handle PDB ID effectively  
def parse_pdb_ids(pdb_ids):
    if isinstance(pdb_ids, list) and len(pdb_ids) > 0:
        return [pdb.strip() for pdb in pdb_ids[0].split(",")]
    return []

# Function to download multiple PDBs
def GetPDB(pdb_ids, folder_path):
    pdbids = parse_pdb_ids(pdb_ids)
    print(pdbids)
    os.makedirs(folder_path, exist_ok=True)
    downloaded_files = []
    for pdb_id in pdbids:
        file_path = os.path.join(folder_path, f"{pdb_id}.pdb")
        try:
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            response = requests.get(url)
            if response.status_code == 200:
                with open(file_path, "w") as f:
                    f.write(response.text)
                downloaded_files.append(file_path)
            else:
                raise Exception(f"Failed to download {pdb_id}.pdb (Status Code: {response.status_code})")
        except Exception as e:
            print(f"Error: {e}")
    return downloaded_files

#Function to download ligand PDB based on PubChem ID
def download_pubchem_sdf(cid_name_pairs, folder_path):
    os.makedirs(folder_path, exist_ok=True)

    for cid, name in cid_name_pairs:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF"
        file_path = os.path.join(folder_path, f"{name}.sdf")
        try:
            response = requests.get(url, stream=True)
            if response.status_code == 200:
                with open(file_path, 'wb') as file:
                    for chunk in response.iter_content(chunk_size=1024):
                        if chunk:
                            file.write(chunk)

                supplier = Chem.SDMolSupplier(file_path)
                mol = next((m for m in supplier if m is not None), None)
                if mol is None:
                    print("Failed to read molecule from SDF.")
                    continue

                AllChem.Compute2DCoords(mol)
                pdb_file = os.path.join(folder_path, f"{name}.pdb")
                Chem.MolToPDBFile(mol, pdb_file)
            else:
                print(f"[✗] Failed for CID {cid}. HTTP Status: {response.status_code}")

        except Exception as e:
            print(f"[!] Error with CID {cid}: {e}")

#Function to convert pdb to pdbqt
def convert_to_pdbqt(input_file, output_file, file_type):
    if file_type == "receptor":
        prepare_script = os.path.join(MGLTOOLS_PATH, "prepare_receptor4.py")
    elif file_type == "ligand":
        prepare_script = os.path.join(MGLTOOLS_PATH, "prepare_ligand4.py")
    else:
        raise ValueError("Invalid file type. Use 'receptor' or 'ligand'.")

    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")

    python_exe = MGL_PYTHON
    if not os.path.exists(python_exe):
        # Fallback to system python2.7
        python_exe = "/usr/bin/python2.7"
        if not os.path.exists(python_exe):
            raise FileNotFoundError(f"Python 2.7 executable not found. Tried: {MGL_PYTHON} and /usr/bin/python2.7")

    command = [
        python_exe,
        prepare_script,
        "-l" if file_type == "ligand" else "-r",
        os.path.abspath(input_file),
        "-o",
        os.path.abspath(output_file)
    ]

    logging.info(f"Running command: {' '.join(command)}")

    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        logging.info(f"Conversion successful. Output: {result.stdout}")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Conversion failed. Error: {e.stderr}")
        return False

#Function to calculate centre for performing Vina Run
def calculate_center(pdbqt_file):
    x_coords, y_coords, z_coords = [], [], []

    with open(pdbqt_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                x_coords.append(float(line[30:38].strip()))
                y_coords.append(float(line[38:46].strip()))
                z_coords.append(float(line[46:54].strip()))

    center_x = round(sum(x_coords) / len(x_coords), 3)
    center_y = round(sum(y_coords) / len(y_coords), 3)
    center_z = round(sum(z_coords) / len(z_coords), 3)

    return center_x, center_y, center_z

#Function to convert pdbqt to pdb
def pdbqt_to_pdb(pdbqt_file, output_file):
    python_exe = MGL_PYTHON
    if not os.path.exists(python_exe):
        python_exe = "/usr/bin/python2.7"
    
    prepare_script = os.path.join(MGLTOOLS_PATH, "pdbqt_to_pdb.py")
    command = [
        python_exe,
        prepare_script,
        "-f",
        os.path.abspath(pdbqt_file),
        "-o",
        os.path.abspath(output_file)
    ]

    logging.info(f"Running command: {' '.join(command)}")

    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        logging.info(f"Conversion successful. Output: {result.stdout}")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Conversion failed. Error: {e.stderr}")
        return False

#Function to make a complex of protein and docked ligand
def protein_ligand_complex_file(protein_pdb, ligand_pdb, output_pdb):
    with open(protein_pdb, 'r') as protein_file, open(ligand_pdb, 'r') as ligand_file, open(output_pdb, 'w') as output_file:
        for line in protein_file:
            if line.startswith("END") or line.startswith("TER"):
                continue  # Ignore END or TER from protein
            elif line.startswith("ATOM"):
                output_file.write(line)

        output_file.write("TER\n")  # Add TER to separate protein and ligand

        for line in ligand_file:
            if line.startswith("END") or line.startswith("TER"):
                continue  # Ignore END or TER from ligand
            elif line.startswith("HETATM"):
                output_file.write(line)

        output_file.write("END\n")  # Final END for the complex

    print(f"Complex PDB file saved as: {output_pdb}")

#Function to make a zip file of all the uploaded and result files
def zip_files(job_id, protein_name, ligand_name, output_file, dirs):
    protein_pdb = os.path.join(dirs['upload_folder'], f"{protein_name}.pdb")
    protein_pdbqt = os.path.join(dirs['protein_pdbqt_dir'], f"{protein_name}.pdbqt")
    ligand_pdb = os.path.join(dirs['upload_folder'], f"{ligand_name}.pdb")
    ligand_pdbqt = os.path.join(dirs['ligand_pdbqt_dir'], f"{ligand_name}.pdbqt")
    docking_pdb = os.path.join(dirs['docking_results_dir'], f"docking_{protein_name}_{ligand_name}.pdb")
    docking_pdbqt = os.path.join(dirs['docking_results_dir'], f"docking_{protein_name}_{ligand_name}.pdbqt")
    complex_pdb = os.path.join(dirs['docking_results_dir'], f"complex_{protein_name}_{ligand_name}.pdb")
    log_file = os.path.join(dirs['docking_results_dir'], f"docking_{protein_name}_{ligand_name}.log")
    
    files_to_zip = [protein_pdb, ligand_pdb, protein_pdbqt, ligand_pdbqt, docking_pdb, docking_pdbqt, complex_pdb, log_file]
    output_zip = output_file
    
    with zipfile.ZipFile(output_zip, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for file in files_to_zip:
            if os.path.exists(file):
                zipf.write(file, arcname=os.path.basename(file))
            else:
                print(f"File not found: {file}")

def extract_energy(log_file_path):
    try:
        with open(log_file_path, 'r') as file:
            for line in file:
                if line.strip().startswith("1"):
                    parts = line.split()
                    if len(parts) >= 2:
                        energy = float(parts[1])
                        return energy
    except Exception as e:
        return None  
    return None

# =========================
# CELERY TASK FOR ASYNC DOCKING
# =========================
@celery.task(bind=True)
def perform_docking_task(self, job_id, protein_files, ligand_files):
    """
    Celery task to perform docking asynchronously
    Now with ISOLATED directories per job
    """
    try:
        # Create isolated directories for this job
        dirs = create_job_directories(job_id)
        
        # Update task state
        self.update_state(state='PROCESSING', meta={'status': 'Converting proteins to PDBQT', 'job_id': job_id})
        
        # Convert proteins to PDBQT
        protein_pdbqts = []
        for protein_path in protein_files:
            protein_pdbqt = os.path.join(dirs['protein_pdbqt_dir'], os.path.basename(protein_path).replace('.pdb', '.pdbqt'))
            if convert_to_pdbqt(protein_path, protein_pdbqt, "receptor"):
                protein_pdbqts.append(protein_pdbqt)

        # Update task state
        self.update_state(state='PROCESSING', meta={'status': 'Converting ligands to PDBQT', 'job_id': job_id})
        
        # Convert ligands to PDBQT
        ligand_pdbqts = []
        for ligand_path in ligand_files:
            ligand_pdbqt = os.path.join(dirs['ligand_pdbqt_dir'], os.path.basename(ligand_path).replace('.pdb', '.pdbqt'))
            if convert_to_pdbqt(ligand_path, ligand_pdbqt, "ligand"):
                ligand_pdbqts.append(ligand_pdbqt)

        # Docking for every protein-ligand pair (ONLY from this job)
        size_x, size_y, size_z = 120, 120, 120
        total_pairs = len(protein_pdbqts) * len(ligand_pdbqts)
        current_pair = 0
        
        results = []
        
        for protein_pdbqt in protein_pdbqts:
            center_x, center_y, center_z = calculate_center(protein_pdbqt)
            for ligand_pdbqt in ligand_pdbqts:
                current_pair += 1
                progress = int((current_pair / total_pairs) * 100)
                
                protein_name = os.path.splitext(os.path.basename(protein_pdbqt))[0]
                ligand_name = os.path.splitext(os.path.basename(ligand_pdbqt))[0]
                
                self.update_state(
                    state='PROCESSING', 
                    meta={
                        'status': f'Docking {protein_name} with {ligand_name}',
                        'progress': progress,
                        'current': current_pair,
                        'total': total_pairs,
                        'job_id': job_id
                    }
                )
                
                output_file_pdbqt = os.path.join(
                    dirs['docking_results_dir'],
                    f"docking_{protein_name}_{ligand_name}.pdbqt"
                )
                log_file = output_file_pdbqt.replace('.pdbqt', '.log')

                command = [
                    VINA_PATH,
                    "--receptor", protein_pdbqt,
                    "--ligand", ligand_pdbqt,
                    "--center_x", str(center_x),
                    "--center_y", str(center_y),
                    "--center_z", str(center_z),
                    "--size_x", str(size_x),
                    "--size_y", str(size_y),
                    "--size_z", str(size_z),
                    "--out", output_file_pdbqt,
                    "--log", log_file
                ]

                try:
                    subprocess.run(command, check=True)
                except subprocess.CalledProcessError as e:
                    logging.error(f'Docking failed for {protein_name}-{ligand_name}: {str(e)}')
                    continue

                output_file_pdb = output_file_pdbqt.replace('.pdbqt', '.pdb')
                pdbqt_to_pdb(output_file_pdbqt, output_file_pdb)

                complex_file = os.path.join(
                    dirs['docking_results_dir'],
                    f"complex_{protein_name}_{ligand_name}.pdb"
                )

                protein_ligand_complex_file(
                    os.path.join(dirs['upload_folder'], f"{protein_name}.pdb"),
                    output_file_pdb, 
                    complex_file
                )

                output_file_zip = output_file_pdbqt.replace('.pdbqt', '.zip')
                zip_files(job_id, protein_name, ligand_name, output_file_zip, dirs)
                
                # Extract binding energy
                energy = extract_energy(log_file)
                
                results.append({
                    'protein': protein_name,
                    'ligand': ligand_name,
                    'energy': energy,
                    'zip_file': os.path.basename(output_file_zip)
                })

        # Cleanup intermediate files but keep results
        cleanup_job_files(job_id, keep_results=True)

        return {
            'status': 'completed',
            'job_id': job_id,
            'results': results,
            'result_folder': dirs['docking_results_dir']
        }
        
    except Exception as e:
        logging.error(f"Docking task failed: {str(e)}")
        # Cleanup on failure
        cleanup_job_files(job_id, keep_results=False)
        return {
            'status': 'failed',
            'job_id': job_id,
            'error': str(e)
        }

# =========================
# FLASK ROUTES
# =========================
@app.route('/')
def index():
    return jsonify({
        'status': 'running',
        'message': 'Molecular Docking Server',
        'version': '2.0'
    })

@app.route('/api/upload', methods=['POST'])
def upload_files():
    """
    Upload files and start async docking job
    NOW WITH ISOLATED DIRECTORIES
    """
    choice = request.form.get('choice')
    protein_files = []
    ligand_files = []
    
    # Generate unique job ID
    job_id = str(uuid.uuid4())
    
    # Create isolated directories for this job
    dirs = create_job_directories(job_id)

    if choice == 'pdb_files':  # Multiple protein & ligand file upload
        proteins = request.files.getlist('proteins')
        ligands = request.files.getlist('ligands')
        if not proteins or not ligands:
            return jsonify({'error': 'Missing protein or ligand files'}), 400

        for protein in proteins:
            # Save to JOB-SPECIFIC upload folder
            protein_path = os.path.join(dirs['upload_folder'], protein.filename)
            protein.save(protein_path)
            protein_files.append(protein_path)

        for ligand in ligands:
            # Save to JOB-SPECIFIC upload folder
            ligand_path = os.path.join(dirs['upload_folder'], ligand.filename)
            ligand.save(ligand_path)
            ligand_files.append(ligand_path)

    elif choice == 'pdb_id_compound':  # Multiple PDB ID input
        pdb_ids = request.form.getlist('pdb_ids')
        ligand_cid_names_raw = request.form.getlist('ligand_cid_name')

        if not pdb_ids or not ligand_cid_names_raw:
            return jsonify({'error': 'Missing PDB IDs or ligand CID-name pairs'}), 400

        # Handle input like ['3314, EUGENOL, 3034034, QUININE']
        ligand_cid_names_flat = []
        for item in ligand_cid_names_raw:
            ligand_cid_names_flat.extend(item.split(','))

        # Now make (CID, Name) pairs
        cid_name_pairs = []
        for i in range(0, len(ligand_cid_names_flat), 2):
            cid = ligand_cid_names_flat[i].strip()
            name = ligand_cid_names_flat[i+1].strip()
            cid_name_pairs.append((cid, name))

        print(cid_name_pairs)

        # Download to JOB-SPECIFIC upload folder
        protein_files = GetPDB(pdb_ids, dirs['upload_folder'])
        download_pubchem_sdf(cid_name_pairs, dirs['upload_folder'])
        ligand_files = [os.path.join(dirs['upload_folder'], f"{name}.pdb") for _, name in cid_name_pairs]

    else:
        return jsonify({'error': 'Invalid choice'}), 400

    # Start async docking task
    task = perform_docking_task.apply_async(args=[job_id, protein_files, ligand_files])
    
    return jsonify({
        'message': 'Docking job started',
        'job_id': job_id,
        'task_id': task.id
    }), 202

@app.route('/api/job_status/<task_id>', methods=['GET'])
def job_status(task_id):
    """
    Check status of docking job
    """
    task = AsyncResult(task_id, app=celery)
    
    if task.state == 'PENDING':
        response = {
            'state': task.state,
            'status': 'Job is waiting in queue'
        }
    elif task.state == 'PROCESSING':
        response = {
            'state': task.state,
            'job_id': task.info.get('job_id', ''),
            'status': task.info.get('status', ''),
            'progress': task.info.get('progress', 0),
            'current': task.info.get('current', 0),
            'total': task.info.get('total', 0)
        }
    elif task.state == 'SUCCESS':
        response = {
            'state': task.state,
            'result': task.result
        }
    else:
        response = {
            'state': task.state,
            'status': str(task.info)
        }
    
    return jsonify(response)

@app.route("/api/download", methods=["GET"])
def download_pdb():
    job_id = request.args.get("job_id")
    protein = request.args.get("protein")
    ligand = request.args.get("ligand")

    if not job_id or not protein or not ligand:
        return jsonify({"error": "Missing job_id, protein, or ligand parameter"}), 400

    file_path = os.path.join(BASE_FOLDER, job_id, "docking_results", f"docking_{protein}_{ligand}.zip")

    if os.path.exists(file_path):
        return send_file(file_path, as_attachment=True)
    else:
        return jsonify({"error": "File not found"}), 404

@app.route("/api/get_complex_file", methods=["GET"])
def get_complex_file():
    job_id = request.args.get("job_id")
    protein = request.args.get("protein")
    ligand = request.args.get("ligand")

    if not job_id or not protein or not ligand:
        return jsonify({"error": "Missing job_id, protein, or ligand parameter"}), 400

    file_path = os.path.join(BASE_FOLDER, job_id, "docking_results", f"complex_{protein}_{ligand}.pdb")
    
    if os.path.exists(file_path):
        return jsonify({"file_name": file_path})
    else:
        return jsonify({"error": "Complex file not found"}), 404
    
@app.route('/api/visualize_pdb', methods=['GET'])
def visualize_pdb():
    file_name = request.args.get("file_name")
    print(f"Received file_name: {file_name}")  
    
    if not file_name:
        print("Error: file_name is None")  
        return jsonify({"error": "file_name parameter is missing"}), 400

    if not os.path.exists(file_name):
        print("Error: File does not exist!")  
        return jsonify({"error": "File not found"}), 404

    with open(file_name, "r") as f:
        pdb_data = f.read()
    
    return jsonify({"pdb_data": pdb_data})

@app.route('/get_energy', methods=['GET'])
def get_energy():
    job_id = request.args.get('job_id')
    protein = request.args.get('protein')
    ligand = request.args.get('ligand')

    if not job_id or not protein or not ligand:
        return jsonify({"error": "Missing job_id, protein, or ligand parameter"}), 400

    log_file_path = os.path.join(BASE_FOLDER, job_id, "docking_results", f"docking_{protein}_{ligand}.log")

    print(f"Looking for file: {log_file_path}")  

    if not os.path.exists(log_file_path):
        print("Log file not found!")  
        return jsonify({"error": "Log file not found"}), 404

    energy = extract_energy(log_file_path)
    print("Energy:", energy)

    if energy is None:
        print("Energy value not found!")  
        return jsonify({"error": "Energy value not found"}), 500

    return jsonify({"interaction_energy": energy})

@app.route('/api/health', methods=['GET'])
def health():
    """Health check endpoint"""
    return jsonify({
        'status': 'healthy',
        'timestamp': datetime.now().isoformat()
    })

@app.route('/api/cleanup_old_jobs', methods=['POST'])
def cleanup_old_jobs():
    """
    Manual cleanup endpoint to remove old job folders
    """
    days_old = request.json.get('days_old', 7)
    
    import time
    current_time = time.time()
    deleted = []
    
    for job_id in os.listdir(BASE_FOLDER):
        job_path = os.path.join(BASE_FOLDER, job_id)
        if os.path.isdir(job_path):
            # Check folder age
            folder_age_days = (current_time - os.path.getmtime(job_path)) / (86400)
            if folder_age_days > days_old:
                shutil.rmtree(job_path)
                deleted.append(job_id)
    
    return jsonify({
        'deleted_jobs': len(deleted),
        'job_ids': deleted
    })

if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5000, debug=False)
