import React from "react";
import { AppBar, Toolbar, Typography, Button } from "@mui/material";
import { motion } from "framer-motion";

const Header = () => {
    return (
        <motion.div initial={{ opacity: 0 }} animate={{ opacity: 1 }} transition={{ duration: 1 }}>
            <AppBar 
                position="fixed" 
                sx={{ 
                    width: "100vw",   
                    background: "linear-gradient(90deg, #0F2027, #203A43, #2C5364)", 
                    top: 0, 
                    left: 0, 
                    zIndex: 1100 
                }}
            >
                <Toolbar sx={{ display: "flex", justifyContent: "space-between" }}>
                    <Typography variant="h6" sx={{ fontWeight: "bold" }}>
                        DockSphere 🌍
                    </Typography>
                </Toolbar>
            </AppBar>
        </motion.div>
    );
};

export default Header;
