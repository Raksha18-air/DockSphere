from flask import Flask, request, jsonify, send_from_directory, render_template, send_file
from flask_cors import CORS
import os
import subprocess
import logging
import requests
import py3Dmol
import zipfile
import glob
from rdkit import Chem
from rdkit.Chem import AllChem


app = Flask(__name__)
CORS(app)

# Paths to required tools and directories
MGLTOOLS_PATH = os.environ.get("MGLTOOLS_PATH", "/opt/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24")
VINA_PATH = os.environ.get("VINA_PATH", "/usr/local/bin/vina")
MGL_PYTHON = os.environ.get("MGL_PYTHON", "/opt/mgltools/bin/python")

UPLOAD_FOLDER = os.environ.get("UPLOAD_FOLDER", "/app/uploads")
RESULT_FOLDER = os.environ.get("RESULT_FOLDER", "/app/results")

PROTEIN_PDBQT_DIR = os.path.join(RESULT_FOLDER, "proteins_pdbqt")
LIGAND_PDBQT_DIR = os.path.join(RESULT_FOLDER, "ligands_pdbqt")
DOCKING_RESULTS_DIR = os.path.join(RESULT_FOLDER, "docking_result")

# Create necessary directories
for directory in [UPLOAD_FOLDER, RESULT_FOLDER, PROTEIN_PDBQT_DIR, LIGAND_PDBQT_DIR, DOCKING_RESULTS_DIR]:
    os.makedirs(directory, exist_ok=True)

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
        raise FileNotFoundError(f"Python executable not found: {python_exe}")

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
def zip_files(protein_name, ligand_name, output_file):
    protein_pdb = os.path.join(UPLOAD_FOLDER, f"{protein_name}.pdb")
    protein_pdbqt = os.path.join(PROTEIN_PDBQT_DIR, f"{protein_name}.pdbqt")
    ligand_pdb = os.path.join(UPLOAD_FOLDER, f"{ligand_name}.pdb")
    ligand_pdbqt = os.path.join(LIGAND_PDBQT_DIR, f"{ligand_name}.pdbqt")
    docking_pdb = os.path.join(DOCKING_RESULTS_DIR, f"docking_{protein_name}_{ligand_name}.pdb")
    docking_pdbqt = os.path.join(DOCKING_RESULTS_DIR, f"docking_{protein_name}_{ligand_name}.pdbqt")
    complex_pdb = os.path.join(DOCKING_RESULTS_DIR, f"complex_{protein_name}_{ligand_name}.pdb")
    log_file = os.path.join(DOCKING_RESULTS_DIR, f"docking_{protein_name}_{ligand_name}.log")
    files_to_zip = [protein_pdb, ligand_pdb, protein_pdbqt, ligand_pdbqt, docking_pdb, docking_pdbqt, complex_pdb, log_file]
    output_zip = output_file
    with zipfile.ZipFile(output_zip, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for file in files_to_zip:
            if os.path.exists(file):
                zipf.write(file, arcname=os.path.basename(file))
            else:
                print(f"File not found: {file}")

#Function to extract energy from log file
def extract_energy(log_file):
    try:
        with open(log_file, 'r') as file:
            for line in file:
                if "-----+------------+----------+----------" in line:
                    next(file)  
                    pose_line = next(file)  
                    energy = float(pose_line.split()[1])  
                    return energy
    except Exception as e:
        return None  
    return None  

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/api/upload', methods=['POST'])
def upload_files():
    choice = request.form.get('choice')
    protein_files = []
    ligand_files = []

    if choice == 'pdb_files':  # Multiple protein & ligand file upload
        proteins = request.files.getlist('proteins')
        ligands = request.files.getlist('ligands')
        if not proteins or not ligands:
            return jsonify({'error': 'Missing protein or ligand files'}), 400

        for protein in proteins:
            protein_path = os.path.join(UPLOAD_FOLDER, protein.filename)
            protein.save(protein_path)
            protein_files.append(protein_path)

        for ligand in ligands:
            ligand_path = os.path.join(UPLOAD_FOLDER, ligand.filename)
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

        protein_files = GetPDB(pdb_ids, UPLOAD_FOLDER)
        download_pubchem_sdf(cid_name_pairs, UPLOAD_FOLDER)
        ligand_files = [os.path.join(UPLOAD_FOLDER, f"{name}.pdb") for _, name in cid_name_pairs]

    else:
        return jsonify({'error': 'Invalid choice'}), 400

    # Convert proteins to PDBQT
    protein_pdbqts = []
    for protein_path in protein_files:
        protein_pdbqt = os.path.join(PROTEIN_PDBQT_DIR, os.path.basename(protein_path).replace('.pdb', '.pdbqt'))
        if convert_to_pdbqt(protein_path, protein_pdbqt, "receptor"):
            protein_pdbqts.append(protein_pdbqt)

    # Convert ligands to PDBQT
    ligand_pdbqts = []
    for ligand_path in ligand_files:
        ligand_pdbqt = os.path.join(LIGAND_PDBQT_DIR, os.path.basename(ligand_path).replace('.pdb', '.pdbqt'))
        if convert_to_pdbqt(ligand_path, ligand_pdbqt, "ligand"):
            ligand_pdbqts.append(ligand_pdbqt)

    # Docking for every protein-ligand pair
    size_x, size_y, size_z = 120, 120, 120
    for protein_pdbqt in protein_pdbqts:
        center_x, center_y, center_z = calculate_center(protein_pdbqt)
        for ligand_pdbqt in ligand_pdbqts:
            protein_name = os.path.splitext(os.path.basename(protein_pdbqt))[0]
            ligand_name = os.path.splitext(os.path.basename(ligand_pdbqt))[0]
            output_file_pdbqt = os.path.join(
                DOCKING_RESULTS_DIR,
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
                return jsonify({'error': f'Docking failed: {str(e)}'}), 500

            output_file_pdb = output_file_pdbqt.replace('.pdbqt', '.pdb')

            pdbqt_to_pdb(output_file_pdbqt, output_file_pdb)

            complex_file = os.path.join(
                DOCKING_RESULTS_DIR,
                f"complex_{protein_name}_{ligand_name}.pdb"
            )

            protein_ligand_complex_file(os.path.join(UPLOAD_FOLDER, f"{protein_name}.pdb"),output_file_pdb, complex_file)

            output_file_zip = output_file_pdbqt.replace('.pdbqt', '.zip')

            zip_files(protein_name, ligand_name, output_file_zip)

    return jsonify({'message': 'Docking completed successfully.', 'files': os.listdir(DOCKING_RESULTS_DIR)})

@app.route("/api/download", methods=["GET"])
def download_pdb():
    protein = request.args.get("protein")
    ligand = request.args.get("ligand")

    if not protein or not ligand:
        return jsonify({"error": "Missing protein or ligand parameter"}), 400

    file_path = os.path.join(DOCKING_RESULTS_DIR, f"docking_{protein}_{ligand}.zip")

    if os.path.exists(file_path):
        return send_file(file_path, as_attachment=True)
    else:
        return jsonify({"error": "File not found"}), 404

@app.route("/api/get_complex_file", methods=["GET"])
def get_complex_file():
    protein = request.args.get("protein")
    ligand = request.args.get("ligand")

    if not protein or not ligand:
        return jsonify({"error": "Missing protein or ligand parameter"}), 400

    complex_file = None
    file_path = os.path.join(DOCKING_RESULTS_DIR, f"complex_{protein}_{ligand}.pdb")
    if os.path.exists(file_path):
        complex_file = file_path
        return jsonify({"file_name": complex_file})
    else:
        return jsonify({"error": "Complex file not found"}), 404
    
@app.route('/api/visualize_pdb', methods=['GET'])
def visualize_pdb():
    file_name = request.args.get("file_name")
    print(f"Received file_name: {file_name}")  
    
    if not file_name:
        print("Error: file_name is None")  
        return jsonify({"error": "file_name parameter is missing"}), 400

    file_path = os.path.join(os.getcwd(), file_name)
    print(f"Constructed file path: {file_path}")  
    
    if not os.path.exists(file_path):
        print("Error: File does not exist!")  
        return jsonify({"error": "File not found"}), 404

    with open(file_path, "r") as f:
        pdb_data = f.read()
    
    return jsonify({"pdb_data": pdb_data})

@app.route('/get_energy', methods=['GET'])
def get_energy():
    protein = request.args.get('protein')
    ligand = request.args.get('ligand')

    log_file_path = os.path.join(os.getcwd(), DOCKING_RESULTS_DIR, f"docking_{protein}_{ligand}.log")

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

@app.route("/get_pdb_data")
def get_pdb_data():
    return jsonify({"pdb": open("results/docking_result/complex_1I78_eugenol.pdb").read()})

if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5000)

