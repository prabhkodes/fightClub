FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install system-level dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    gfortran \
    valgrind \
    gdb \
    cmake \
    vim \
    mpich \
    libpnetcdf-dev \
    libnetcdf-dev \
    libnetcdff-dev \
    python3 \
    python3-pip \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install Python dependencies
COPY requirements.txt /tmp/requirements.txt
RUN pip3 install --no-cache-dir -r /tmp/requirements.txt

# Workspace
WORKDIR /app

CMD ["/bin/bash"]
