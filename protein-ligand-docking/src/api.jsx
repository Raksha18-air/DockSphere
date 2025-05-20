import axios from "axios";

const API_BASE_URL = "http://127.0.0.1:5000/api"; // Adjust if needed

export const uploadFiles = (formData) => {
    return axios.post(`${API_BASE_URL}/upload`, formData);
};

export const fetchResults = () => {
    return axios.get(`${API_BASE_URL}/results`);
};

export const fetchVisualization = () => {
    return axios.get(`${API_BASE_URL}/visualize`);
};
