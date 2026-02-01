import React, { useState, useContext } from "react";
import { Box, Button, Select, MenuItem, TextField, Typography, List, ListItem, ListItemText, Container } from "@mui/material";
import { motion } from "framer-motion";
import { DockingContext } from "../DockingContext";  

const UploadForm = () => {
    const { setDockingResults, setIsLoading } = useContext(DockingContext);
    const [inputType, setInputType] = useState("pdb_files");
    const [pdbId, setPdbId] = useState("");
    const [ligandCid, setLigandCid] = useState("");
    const [proteinFiles, setProteinFiles] = useState([]);
    const [ligandFiles, setLigandFiles] = useState([]);

    const handleProteinFileChange = (event) => {
        setProteinFiles(Array.from(event.target.files));
    };

    const handleLigandFileChange = (event) => {
        setLigandFiles(Array.from(event.target.files));
    };

    const isSubmitDisabled = inputType === "pdb_files" && (proteinFiles.length === 0 || ligandFiles.length === 0);

    const API_BASE_URL = process.env.REACT_APP_API_URL || "http://127.0.0.1:5000";

    const handleSubmit = async () => {
        setIsLoading(true);
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
            const response = await fetch(`${API_BASE_URL}/api/upload`, {
                method: "POST",
                body: formData
            });

            const result = await response.json();
            setIsLoading(false);
            setDockingResults(result);  

            if (response.ok) {
                alert("Docking Successful!");
            } else {
                alert("Error: " + result.error);
            }
        } catch (error) {
            console.error("Error uploading files:", error);
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
                            <Typography variant="subtitle1">Upload Protein PDB Files (Don't use _ in file name and upload PDB without hetroatoms):</Typography>
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
                            <TextField label="Ligand CID" placeholder="Enter Ligand CID(s) and Name(s) (Comma Seperated) (e.g., 71853, Chitosan)" fullWidth margin="normal" value={ligandCid} onChange={(e) => setLigandCid(e.target.value)} />
                        </>
                    )}

                    <Button fullWidth variant="contained" sx={{ mt: 2, background: "linear-gradient(90deg, #ff416c, #ff4b2b)" }}
                        disabled={isSubmitDisabled} onClick={handleSubmit}>
                        Submit & Start Docking
                    </Button>
                </Container>
            </Box>
        </motion.div>
    );
};

export default UploadForm;
