import React, { useContext } from "react";
import axios from "axios";
import { DockingContext } from "../DockingContext";
import { Table, TableBody, TableCell, TableContainer, TableHead, TableRow, Button, Paper, Typography  } from "@mui/material";

const ResultsTable = () => {
    const { dockingResults, setSelectedFile } = useContext(DockingContext);

    const results = dockingResults || { files: [] };


    console.log("Docking Results from Context:", dockingResults);

    const API_BASE_URL = process.env.REACT_APP_API_URL || "http://127.0.0.1:5000";

    const handleDownload = (protein, ligand) => {
        console.log(`Downloading file for: Protein - ${protein}, Ligand - ${ligand}`);
        
        axios.get(`${API_BASE_URL}/download?protein=${protein}&ligand=${ligand}`, { responseType: "blob" })
            .then(response => {
                console.log("File downloaded successfully!");
                
                const url = window.URL.createObjectURL(new Blob([response.data]));
                const link = document.createElement("a");
                link.href = url;
                link.setAttribute("download", `${protein}_${ligand}.zip`);
                document.body.appendChild(link);
                link.click();
                link.remove();
            })
            .catch(error => console.error("Download error:", error));
    };
    
    const handleVisualize = (protein, ligand) => {
        console.log(`Fetching complex file for: Protein - ${protein}, Ligand - ${ligand}`);
        
        axios.get(`${API_BASE_URL}/get_complex_file?protein=${protein}&ligand=${ligand}`)
            .then(response => {
                console.log("Received response for visualization:", response.data);

                if (response.data.file_name) {
                    setSelectedFile(response.data.file_name);
                    console.log("Selected file set for visualization:", response.data.file_name);
                } else {
                    console.error("Complex file not found!");
                }
            })
            .catch(error => console.error("Error fetching file:", error));
    };


    if (!dockingResults) {
        return <Typography variant="h6" align="center" color="textSecondary" style={{ marginTop: "20px" }}>
            Stay tuned! Exciting docking results and interactive visualizations are on their way...
        </Typography>;
    }

    return (
        <div>
            <TableContainer component={Paper}>
                <Table>
                    <TableHead>
                        <TableRow>  
                            <TableCell>Protein</TableCell>
                            <TableCell>Ligand</TableCell>
                            <TableCell>Actions</TableCell>
                        </TableRow>
                    </TableHead>
                    <TableBody>
                        {Array.isArray(dockingResults.files) && dockingResults.files.length > 0 ? (
                            (() => {
                                const complexFiles = dockingResults.files.filter(file => file.startsWith("complex"));

                                if (complexFiles.length === 0) {
                                    return (
                                        <TableRow>
                                            <TableCell colSpan={3} style={{ textAlign: "center", fontWeight: "bold", color: "red" }}>
                                                No valid docking results found.
                                            </TableCell>
                                        </TableRow>
                                    );
                                }

                                return complexFiles.map((file, index) => {
                                    const fileParts = file.replace(".pdb", "").split("_");

                                    if (fileParts.length >= 3) {
                                        const protein = fileParts[1];
                                        const ligand = fileParts[2];

                                        return (
                                            <TableRow key={index}>
                                                <TableCell>{protein}</TableCell>
                                                <TableCell>{ligand}</TableCell>
                                                <TableCell>
                                                    <Button variant="contained" color="primary" onClick={() => handleDownload(protein, ligand)}>
                                                        Download
                                                    </Button>
                                                    &nbsp;
                                                    <Button variant="contained" color="secondary" onClick={() => handleVisualize(protein, ligand)}>
                                                        Visualize
                                                    </Button>
                                                </TableCell>
                                            </TableRow>
                                        );
                                    } else {
                                        return (
                                            <TableRow key={index}>
                                                <TableCell colSpan={3} style={{ textAlign: "center", color: "red" }}>
                                                    Invalid file format: {file}
                                                </TableCell>
                                            </TableRow>
                                        );
                                    }
                                });
                            })()
                        ) : (
                            <TableRow>
                                <TableCell colSpan={3} style={{ textAlign: "center", fontWeight: "bold", color: "red" }}>
                                    No docking results found.
                                </TableCell>
                            </TableRow>
                        )}
                    </TableBody>
                </Table>
            </TableContainer>

        </div>
    );
};

export default ResultsTable;
