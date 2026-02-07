FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# =========================
# System Dependencies
# =========================
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python2.7 \
    wget \
    curl \
    unzip \
    build-essential \
    libgl1 \
    libglu1-mesa \
    libxrender1 \
    libxext6 \
    libxi6 \
    ca-certificates \
    redis-server \
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
# Download MGLTools directly if not provided locally
#RUN wget -q https://ccsb.scripps.edu/mgltools/download/491/ -O /tmp/mgltools.tar.gz || \
#echo "Please ensure mgltools_x86_64Linux2_1.5.7.tar.gz is available"

# If you have the file locally, use this instead:
COPY backend/mgltools_x86_64Linux2_*.tar.gz /tmp/mgltools.tar.gz

RUN mkdir -p /opt && \
    tar -xzf /tmp/mgltools.tar.gz -C /opt/ && \
    mv /opt/mgltools_x86_64Linux2_* /opt/mgltools && \
    rm /tmp/mgltools.tar.gz

# Run the MGLTools installation script
RUN cd /opt/mgltools && \
    ./install.sh -d /opt/mgltools -c 1

# Verify Python 2.7 is available in MGLTools
RUN ls -la /opt/mgltools/bin/ && \
    test -f /opt/mgltools/bin/python2.7 || \
    ln -s /usr/bin/python2.7 /opt/mgltools/bin/python

ENV MGLTOOLS_PATH=/opt/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24
ENV MGL_PYTHON=/opt/mgltools/bin/python2.7
ENV VINA_PATH=/usr/local/bin/vina
ENV PATH="/opt/mgltools/bin:${PATH}"

# =========================
# Backend App Setup
# =========================
WORKDIR /app

COPY backend/requirements.txt /app/requirements.txt
RUN pip3 install --no-cache-dir -r /app/requirements.txt

# Install Celery for async job processing
RUN pip3 install celery redis flower

COPY backend/ /app/

# Create necessary directories
RUN mkdir -p /app/uploads /app/results /app/results/proteins_pdbqt \
    /app/results/ligands_pdbqt /app/results/docking_result

# =========================
# Expose Ports
# =========================
EXPOSE 5000 5555

# =========================
# Health Check
# =========================
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:5000/ || exit 1

# Start script will run both Flask and Celery
COPY backend/start.sh /app/start.sh
RUN chmod +x /app/start.sh

CMD ["/app/start.sh"]
