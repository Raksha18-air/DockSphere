import React from "react";
import { AppBar, Toolbar, Typography, IconButton } from "@mui/material";
import { motion } from "framer-motion";
import GitHubIcon from "@mui/icons-material/GitHub";
import EmailIcon from "@mui/icons-material/Email";
import LinkedInIcon from "@mui/icons-material/LinkedIn";

const Footer = () => {
    return (
        <motion.div initial={{ opacity: 0 }} animate={{ opacity: 1 }} transition={{ duration: 1 }}>
            <AppBar 
                position="static" 
                sx={{ 
                    background: "linear-gradient(90deg, #0F2027, #203A43, #2C5364)", 
                    top: "auto", 
                    bottom: 0, 
                    marginTop: "2rem"
                }}
            >
                <Toolbar sx={{ flexDirection: "column", textAlign: "center", py: 2 }}>
                    <Typography variant="body2" color="inherit">
                        &copy; 2025 DockSphere 🌍 | Managed by Raksha Ray.
                    </Typography>
                    <Typography variant="body2" color="inherit">
                        Designed and Developed by Raksha Ray
                    </Typography>

                    <div>
                        <IconButton href="https://github.com/Raksha18-air" target="_blank" color="inherit">
                            <GitHubIcon />
                        </IconButton>
                        <IconButton href="mailto:r.raksha@iitg.ac.in" color="inherit">
                            <EmailIcon />
                        </IconButton>
                        <IconButton href="https://www.linkedin.com/in/raksha-ray-1055271a9/" target="_blank" color="inherit">
                            <LinkedInIcon />
                        </IconButton>
                    </div>
                </Toolbar>
            </AppBar>
        </motion.div>
    );
};

export default Footer;
