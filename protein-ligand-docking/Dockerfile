# =========================
# Build Stage
# =========================
FROM node:18-alpine AS build

WORKDIR /frontend

# Copy env for build-time
COPY protein-ligand-docking/.env .env

COPY protein-ligand-docking/package*.json ./
RUN npm install

COPY protein-ligand-docking/ ./
RUN npm run build

# =========================
# Production Stage (Nginx)
# =========================
FROM nginx:alpine

# Remove default nginx static files
RUN rm -rf /usr/share/nginx/html/*

# Copy React build to Nginx
COPY --from=build /frontend/build /usr/share/nginx/html

EXPOSE 80

CMD ["nginx", "-g", "daemon off;"]