import React, { useEffect, useState, useContext, useRef } from "react";
import axios from "axios";
import { Card, CardContent, Typography, CircularProgress, IconButton } from "@mui/material";
import ZoomInIcon from "@mui/icons-material/ZoomIn";
import ZoomOutIcon from "@mui/icons-material/ZoomOut";
import RefreshIcon from "@mui/icons-material/Refresh";
import FullscreenIcon from "@mui/icons-material/Fullscreen";
import { DockingContext } from "../DockingContext";

const Visualization = () => {
    const { selectedFile } = useContext(DockingContext);
    const [loading, setLoading] = useState(false);
    const [pdbData, setPdbData] = useState("");
    const viewerContainerRef = useRef(null);
    const viewerRef = useRef(null);

    useEffect(() => {
        if (selectedFile) {
            console.log("Selected file:", selectedFile);
            setLoading(true);

            const encodedFileName = encodeURIComponent(selectedFile);

            axios.get(`http://127.0.0.1:5000/api/visualize_pdb?file_name=${encodedFileName}`)
                .then(response => {
                    console.log("API Response:", response);
                    if (response.data.pdb_data) {
                        console.log("PDB Data received!", response.data.pdb_data);
                        setPdbData(response.data.pdb_data);
                    } else {
                        console.error("PDB Data missing in response:", response.data);
                    }
                })
                .catch(error => {
                    console.error("Error fetching PDB data:", error);
                })
                .finally(() => {
                    setLoading(false);
                    console.log("Loading state set to false");
                });
        }
    }, [selectedFile]);

    useEffect(() => {
        if (pdbData && viewerContainerRef.current) {
            console.log("Initializing 3Dmol.js viewer...");

            try {
                const container = viewerContainerRef.current;
                if (!container) {
                    console.error("Viewer container not found!");
                    return;
                }

                // Purana viewer clear karo agar exist karta hai
                if (viewerRef.current) {
                    viewerRef.current.clear();
                    viewerRef.current = null;
                }
                viewerRef.current = window.$3Dmol.createViewer(container, { backgroundColor: "white" });

                let viewer = viewerRef.current;
                viewer.clear();

                let model = viewer.addModel(pdbData, "pdb");
                console.log("Applying visualization style...");
                model.setStyle({}, { line: { color: "green", linewidth: 1.5 } });
                model.setStyle({ hetflag: true }, { stick: { radius: 0.3 } });

                viewer.zoomTo();
                viewer.render();
                console.log("Viewer render complete!");
            } catch (error) {
                console.error("Error in 3Dmol.js visualization:", error);
            }
        }
    }, [pdbData]);

    // 🔍 Zoom In
    const zoomIn = () => {
        if (viewerRef.current) {
            viewerRef.current.zoom(1.2);
            viewerRef.current.render();
        }
    };

    // 🔍 Zoom Out
    const zoomOut = () => {
        if (viewerRef.current) {
            viewerRef.current.zoom(0.8);
            viewerRef.current.render();
        }
    };

    // 🔄 Reset View
    const resetView = () => {
        if (viewerRef.current) {
            viewerRef.current.zoomTo();
            viewerRef.current.render();
        }
    };

    // 🖥 Fullscreen Toggle
    const toggleFullscreen = () => {
        const elem = viewerContainerRef.current;
        if (elem) {
            if (!document.fullscreenElement) {
                elem.requestFullscreen();
            } else {
                document.exitFullscreen();
            }
        }
    };

    return (
        <Card sx={{ maxWidth: 800, margin: "auto", mt: 4, p: 2, textAlign: "center", position: "relative" }}>
            <CardContent>
                <Typography variant="h5" gutterBottom>
                    Protein-Ligand 3D Visualization
                </Typography>

                {loading ? (
                    <CircularProgress />
                ) : selectedFile && pdbData ? (
                    <div 
                        style={{
                            position: "relative",
                            display: "flex",
                            justifyContent: "center"
                        }}
                    >
                        <div
                            ref={viewerContainerRef}
                            id="viewerContainer"
                            style={{
                                height: "500px",
                                width: "100%",
                                border: "1px solid #ddd",
                                borderRadius: "8px",
                                overflow: "hidden",
                                position: "relative"
                            }}
                        ></div>

                        {/* 🔘 Small Floating Buttons */}
                        <div
                            style={{
                                position: "absolute",
                                top: "10px",
                                right: "10px",
                                display: "flex",
                                flexDirection: "column",
                                background: "rgba(255, 255, 255, 0.8)",
                                borderRadius: "8px",
                                padding: "5px"
                            }}
                        >
                            <IconButton onClick={zoomIn} size="small" title="Zoom In">
                                <ZoomInIcon />
                            </IconButton>
                            <IconButton onClick={zoomOut} size="small" title="Zoom Out">
                                <ZoomOutIcon />
                            </IconButton>
                            <IconButton onClick={resetView} size="small" title="Reset View">
                                <RefreshIcon />
                            </IconButton>
                            <IconButton onClick={toggleFullscreen} size="small" title="Fullscreen">
                                <FullscreenIcon />
                            </IconButton>
                        </div>
                    </div>
                ) : (
                    <Typography>No file selected for visualization.</Typography>
                )}
            </CardContent>
        </Card>
    );
};

export default Visualization;
