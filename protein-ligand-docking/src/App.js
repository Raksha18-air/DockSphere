import React, { useState } from "react";
import UploadForm from "./components/UploadForm";
import ResultsTable from "./components/ResultsTable";
import Visualization from "./components/Visualization";
import { Container, Box } from "@mui/material";
import Header from "./components/Header";
import Footer from "./components/Footer";
import { DockingProvider } from "./DockingContext";

function App() {
    return (
        <DockingProvider>
            <Header />
            <Container>
                <Box sx={{ marginTop: 10 }}>
                    <UploadForm />
                    <ResultsTable />
                    <Visualization />
                    </Box>
                </Container>
            <Footer />
        </DockingProvider>
    );
}

export default App;

