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

COPY protein-ligand-docking/nginx.conf /etc/nginx/conf.d/default.conf
COPY --from=build /app/build /usr/share/nginx/html

EXPOSE 80

CMD ["nginx", "-g", "daemon off;"]
