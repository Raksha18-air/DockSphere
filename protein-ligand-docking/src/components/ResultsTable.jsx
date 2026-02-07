import React, { useContext } from "react";
import axios from "axios";
import { DockingContext } from "../DockingContext";
import { Table, TableBody, TableCell, TableContainer, TableHead, TableRow, Button, Paper, Typography, Chip  } from "@mui/material";

const ResultsTable = () => {
    const { dockingResults, setSelectedFile } = useContext(DockingContext);

    console.log("Docking Results from Context:", dockingResults);

    const API_BASE_URL = process.env.REACT_APP_API_URL || "http://127.0.0.1:5000";

    // Extract job_id from docking results
    const jobId = dockingResults?.job_id;
    const results = dockingResults?.results || [];

    const handleDownload = (protein, ligand) => {
        if (!jobId) {
            alert("Job ID not found! Cannot download results.");
            return;
        }

        console.log(`Downloading file for: Job - ${jobId}, Protein - ${protein}, Ligand - ${ligand}`);
        
        // Updated API call with job_id
        axios.get(`${API_BASE_URL}/download?job_id=${jobId}&protein=${protein}&ligand=${ligand}`, { 
            responseType: "blob" 
        })
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
        .catch(error => {
            console.error("Download error:", error);
            alert("Failed to download file. Please try again.");
        });
    };
    
    const handleVisualize = (protein, ligand) => {
        if (!jobId) {
            alert("Job ID not found! Cannot visualize results.");
            return;
        }

        console.log(`Fetching complex file for: Job - ${jobId}, Protein - ${protein}, Ligand - ${ligand}`);
        
        // Updated API call with job_id
        axios.get(`${API_BASE_URL}/get_complex_file?job_id=${jobId}&protein=${protein}&ligand=${ligand}`)
            .then(response => {
                console.log("Received response for visualization:", response.data);

                if (response.data.file_name) {
                    setSelectedFile(response.data.file_name);
                    console.log("Selected file set for visualization:", response.data.file_name);
                } else {
                    console.error("Complex file not found!");
                    alert("Visualization file not found!");
                }
            })
            .catch(error => {
                console.error("Error fetching file:", error);
                alert("Failed to fetch visualization file. Please try again.");
            });
    };

    if (!dockingResults) {
        return <Typography variant="h6" align="center" color="textSecondary" style={{ marginTop: "20px" }}>
            Stay tuned! Exciting docking results and interactive visualizations are on their way...
        </Typography>;
    }

    // Handle legacy format (files array) and new format (results array)
    const displayResults = () => {
        // New format: results array with protein, ligand, energy
        if (Array.isArray(results) && results.length > 0) {
            return results.map((result, index) => (
                <TableRow key={index}>
                    <TableCell>{result.protein}</TableCell>
                    <TableCell>{result.ligand}</TableCell>
                    <TableCell>
                        {result.energy ? (
                            <Chip 
                                label={`${result.energy} kcal/mol`} 
                                color={result.energy < -7 ? "success" : result.energy < -5 ? "primary" : "default"}
                                size="small"
                            />
                        ) : (
                            <Typography variant="body2" color="textSecondary">N/A</Typography>
                        )}
                    </TableCell>
                    <TableCell>
                        <Button 
                            variant="contained" 
                            color="primary" 
                            size="small"
                            onClick={() => handleDownload(result.protein, result.ligand)}
                            sx={{ mr: 1 }}
                        >
                            Download
                        </Button>
                        <Button 
                            variant="contained" 
                            color="secondary" 
                            size="small"
                            onClick={() => handleVisualize(result.protein, result.ligand)}
                        >
                            Visualize
                        </Button>
                    </TableCell>
                </TableRow>
            ));
        }
        
        // Legacy format: files array (for backward compatibility)
        if (Array.isArray(dockingResults.files) && dockingResults.files.length > 0) {
            const complexFiles = dockingResults.files.filter(file => file.startsWith("complex"));

            if (complexFiles.length === 0) {
                return (
                    <TableRow>
                        <TableCell colSpan={4} style={{ textAlign: "center", fontWeight: "bold", color: "red" }}>
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
                                <Typography variant="body2" color="textSecondary">N/A</Typography>
                            </TableCell>
                            <TableCell>
                                <Button 
                                    variant="contained" 
                                    color="primary" 
                                    size="small"
                                    onClick={() => handleDownload(protein, ligand)}
                                    sx={{ mr: 1 }}
                                >
                                    Download
                                </Button>
                                <Button 
                                    variant="contained" 
                                    color="secondary" 
                                    size="small"
                                    onClick={() => handleVisualize(protein, ligand)}
                                >
                                    Visualize
                                </Button>
                            </TableCell>
                        </TableRow>
                    );
                } else {
                    return (
                        <TableRow key={index}>
                            <TableCell colSpan={4} style={{ textAlign: "center", color: "red" }}>
                                Invalid file format: {file}
                            </TableCell>
                        </TableRow>
                    );
                }
            });
        }

        // No results
        return (
            <TableRow>
                <TableCell colSpan={4} style={{ textAlign: "center", fontWeight: "bold", color: "red" }}>
                    No docking results found.
                </TableCell>
            </TableRow>
        );
    };

    return (
        <div>
            {/* Display Job ID */}
            {jobId && (
                <Paper sx={{ p: 2, mb: 2, backgroundColor: "#f5f5f5" }}>
                    <Typography variant="body2" color="textSecondary">
                        <strong>Job ID:</strong> {jobId}
                    </Typography>
                    <Typography variant="caption" color="textSecondary" display="block">
                        Use this ID to retrieve your results later
                    </Typography>
                </Paper>
            )}

            <TableContainer component={Paper}>
                <Table>
                    <TableHead>
                        <TableRow>  
                            <TableCell><strong>Protein</strong></TableCell>
                            <TableCell><strong>Ligand</strong></TableCell>
                            <TableCell><strong>Binding Energy</strong></TableCell>
                            <TableCell><strong>Actions</strong></TableCell>
                        </TableRow>
                    </TableHead>
                    <TableBody>
                        {displayResults()}
                    </TableBody>
                </Table>
            </TableContainer>
        </div>
    );
};

export default ResultsTable;
