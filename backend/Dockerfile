FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# =========================
# System Dependencies
# =========================
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    wget \
    unzip \
    build-essential \
    libgl1 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# =========================
# Install AutoDock Vina
# =========================
RUN wget -q https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 \
    -O /usr/local/bin/vina \
    && chmod +x /usr/local/bin/vina

# =========================
# Install MGLTools 1.5.7
# =========================
RUN wget -q http://mgltools.scripps.edu/downloads/tars/releases/REL1.5.7/mgltools_x86_64Linux2_1.5.7.tar.gz \
    && tar -xzf mgltools_x86_64Linux2_1.5.7.tar.gz \
    && mv mgltools_x86_64Linux2_1.5.7 /opt/mgltools

ENV MGLTOOLS_PATH=/opt/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24
ENV MGL_PYTHON=/opt/mgltools/bin/python
ENV VINA_PATH=/usr/local/bin/vina

# =========================
# Backend App Setup
# =========================
WORKDIR /app

COPY backend/requirements.txt /app/requirements.txt
RUN pip3 install --no-cache-dir -r /app/requirements.txt

COPY backend/ /app/

# =========================
# Runtime Folders
# =========================
RUN mkdir -p /app/uploads /app/results

EXPOSE 5000

CMD ["python3", "app.py"]