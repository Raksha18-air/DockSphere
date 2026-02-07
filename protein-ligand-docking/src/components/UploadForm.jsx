import React, { useState, useContext } from "react";
import { Box, Button, Select, MenuItem, TextField, Typography, List, ListItem, ListItemText, Container, LinearProgress, Alert } from "@mui/material";
import { motion } from "framer-motion";
import { DockingContext } from "../DockingContext";  

const UploadForm = () => {
    const { setDockingResults, setIsLoading } = useContext(DockingContext);
    const [inputType, setInputType] = useState("pdb_files");
    const [pdbId, setPdbId] = useState("");
    const [ligandCid, setLigandCid] = useState("");
    const [proteinFiles, setProteinFiles] = useState([]);
    const [ligandFiles, setLigandFiles] = useState([]);
    
    // New state for async job tracking
    const [jobStatus, setJobStatus] = useState(null); // 'PENDING', 'PROCESSING', 'SUCCESS', 'FAILURE'
    const [jobProgress, setJobProgress] = useState(0);
    const [currentTaskId, setCurrentTaskId] = useState(null);
    const [currentJobId, setCurrentJobId] = useState(null);
    const [statusMessage, setStatusMessage] = useState("");

    const handleProteinFileChange = (event) => {
        setProteinFiles(Array.from(event.target.files));
    };

    const handleLigandFileChange = (event) => {
        setLigandFiles(Array.from(event.target.files));
    };

    const isSubmitDisabled = inputType === "pdb_files" && (proteinFiles.length === 0 || ligandFiles.length === 0);

    const API_BASE_URL = process.env.REACT_APP_API_URL || "http://127.0.0.1:5000";

    // Poll for job status
    const pollJobStatus = (taskId, jobId) => {
        const pollInterval = setInterval(async () => {
            try {
                const response = await fetch(`${API_BASE_URL}/job_status/${taskId}`);
                const status = await response.json();

                console.log("Job Status:", status);

                if (status.state === 'PENDING') {
                    setJobStatus('PENDING');
                    setStatusMessage('Job is waiting in queue...');
                    setJobProgress(0);
                } else if (status.state === 'PROCESSING') {
                    setJobStatus('PROCESSING');
                    setStatusMessage(status.status || 'Processing...');
                    setJobProgress(status.progress || 0);
                } else if (status.state === 'SUCCESS') {
                    // Job completed successfully
                    clearInterval(pollInterval);
                    setJobStatus('SUCCESS');
                    setStatusMessage('Docking completed successfully!');
                    setJobProgress(100);
                    setIsLoading(false);

                    // Store job_id for later use (download/visualization)
                    const resultWithJobId = {
                        ...status.result,
                        job_id: jobId
                    };
                    
                    setDockingResults(resultWithJobId);
                    alert("Docking Successful!");
                } else if (status.state === 'FAILURE') {
                    // Job failed
                    clearInterval(pollInterval);
                    setJobStatus('FAILURE');
                    setStatusMessage('Docking failed. Please try again.');
                    setIsLoading(false);
                    alert("Error: " + (status.status || "Unknown error occurred"));
                }
            } catch (error) {
                console.error("Error polling job status:", error);
                clearInterval(pollInterval);
                setJobStatus('FAILURE');
                setStatusMessage('Error checking job status');
                setIsLoading(false);
            }
        }, 2000); // Poll every 2 seconds

        // Store interval ID to clear it later if needed
        return pollInterval;
    };

    const handleSubmit = async () => {
        setIsLoading(true);
        setJobStatus('PENDING');
        setStatusMessage('Submitting job...');
        setJobProgress(0);

        const formData = new FormData();
        formData.append("choice", inputType);

        if (inputType === "pdb_files") {
            proteinFiles.forEach(file => formData.append("proteins", file));  
            ligandFiles.forEach(file => formData.append("ligands", file));  
        } else if (inputType === "pdb_id_compound") {
            formData.append("pdb_ids", pdbId);
            formData.append("ligand_cid_name", ligandCid);
        }

        try {
            // Submit job (now returns immediately with job_id and task_id)
            const response = await fetch(`${API_BASE_URL}/upload`, {
                method: "POST",
                body: formData
            });

            const result = await response.json();
            console.log("Job submitted:", result);

            if (response.ok && result.task_id && result.job_id) {
                // Store task_id and job_id
                setCurrentTaskId(result.task_id);
                setCurrentJobId(result.job_id);
                
                // Start polling for status
                pollJobStatus(result.task_id, result.job_id);
                
                setStatusMessage('Job submitted successfully! Processing...');
            } else {
                setIsLoading(false);
                setJobStatus('FAILURE');
                alert("Error: " + (result.error || "Failed to submit job"));
            }
        } catch (error) {
            console.error("Error uploading files:", error);
            setIsLoading(false);
            setJobStatus('FAILURE');
            setStatusMessage('Failed to upload. Please try again!');
            alert("Failed to upload. Please try again!");
        }
    };

    return (
        <motion.div initial={{ y: -20, opacity: 0 }} animate={{ y: 0, opacity: 1 }} transition={{ duration: 1 }}>
            <Box sx={{ display: "flex", justifyContent: "center", alignItems: "center", minHeight: "90vh", mt: -4 }}>
                <Container sx={{ width: "100%", maxWidth: "500px", bgcolor: "white", boxShadow: 3, borderRadius: 2, p: 4, textAlign: "center" }}>
                    <Typography variant="h3" align="center" gutterBottom>
                        DockSphere 🌍
                    </Typography>

                    <Typography variant="body1" align="center" gutterBottom>
                        This docking app, based on Autodock Vina, allows you to upload protein and ligand PDB files or enter PDB IDs and PubChem CIDs. You can provide multiple entries, and the app will perform docking and provide results with 3D visualization.
                    </Typography>
                    <Typography variant="h5" align="center" gutterBottom>
                         Choose an option
                    </Typography>

                    <Select fullWidth value={inputType} onChange={(e) => setInputType(e.target.value)} sx={{ mb: 2 }}>
                        <MenuItem value="pdb_files">Upload PDB Files </MenuItem>
                        <MenuItem value="pdb_id_compound">Enter PDB IDs and CIDs</MenuItem>
                    </Select>

                    {inputType === "pdb_files" && (
                        <>
                            <Typography variant="subtitle1">Upload Protein PDB Files (Don't use _ in file name and upload PDB without heteroatoms):</Typography>
                            <input type="file" accept=".pdb" multiple onChange={handleProteinFileChange} style={{ marginBottom: "10px", width: "100%" }} />
                            <List>{proteinFiles.map((file, index) => <ListItem key={index}><ListItemText primary={file.name} /></ListItem>)}</List>

                            <Typography variant="subtitle1">Upload Ligand PDB Files (Don't use _ in file name):</Typography>
                            <input type="file" accept=".pdb" multiple onChange={handleLigandFileChange} style={{ marginBottom: "10px", width: "100%" }} />
                            <List>{ligandFiles.map((file, index) => <ListItem key={index}><ListItemText primary={file.name} /></ListItem>)}</List>
                        </>
                    )}

                    {inputType === "pdb_id_compound" && (
                        <>
                            <TextField label="PDB ID" placeholder="Enter PDB ID(s) (e.g., 6LZG)" fullWidth margin="normal" value={pdbId} onChange={(e) => setPdbId(e.target.value)} />
                            <TextField label="Ligand CID" placeholder="Enter Ligand CID(s) and Name(s) (Comma Separated) (e.g., 71853, Chitosan)" fullWidth margin="normal" value={ligandCid} onChange={(e) => setLigandCid(e.target.value)} />
                        </>
                    )}

                    <Button fullWidth variant="contained" sx={{ mt: 2, background: "linear-gradient(90deg, #ff416c, #ff4b2b)" }}
                        disabled={isSubmitDisabled || jobStatus === 'PROCESSING'} onClick={handleSubmit}>
                        {jobStatus === 'PROCESSING' ? 'Processing...' : 'Submit & Start Docking'}
                    </Button>

                    {/* Job Status Display */}
                    {jobStatus && (
                        <Box sx={{ mt: 3 }}>
                            {jobStatus === 'PENDING' && (
                                <Alert severity="info">
                                    <Typography variant="body2">{statusMessage}</Typography>
                                </Alert>
                            )}
                            
                            {jobStatus === 'PROCESSING' && (
                                <>
                                    <Alert severity="info">
                                        <Typography variant="body2">{statusMessage}</Typography>
                                    </Alert>
                                    <Box sx={{ mt: 2 }}>
                                        <LinearProgress variant="determinate" value={jobProgress} />
                                        <Typography variant="caption" sx={{ mt: 1 }}>
                                            Progress: {jobProgress}%
                                        </Typography>
                                    </Box>
                                </>
                            )}
                            
                            {jobStatus === 'SUCCESS' && (
                                <Alert severity="success">
                                    <Typography variant="body2">{statusMessage}</Typography>
                                    <Typography variant="caption" display="block" sx={{ mt: 1 }}>
                                        Job ID: {currentJobId}
                                    </Typography>
                                </Alert>
                            )}
                            
                            {jobStatus === 'FAILURE' && (
                                <Alert severity="error">
                                    <Typography variant="body2">{statusMessage}</Typography>
                                </Alert>
                            )}
                        </Box>
                    )}
                </Container>
            </Box>
        </motion.div>
    );
};

export default UploadForm;
