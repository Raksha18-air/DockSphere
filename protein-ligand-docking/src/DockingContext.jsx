import { createContext, useState } from "react";

export const DockingContext = createContext();

export const DockingProvider = ({ children }) => {
    const [dockingResults, setDockingResults] = useState(null);
    const [isLoading, setIsLoading] = useState(false);
    const [selectedFile, setSelectedFile] = useState(null);  

    return (
        <DockingContext.Provider value={{ dockingResults, setDockingResults, isLoading, setIsLoading, selectedFile, setSelectedFile }}>
            {children}
        </DockingContext.Provider>
    );
};
