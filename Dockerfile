FROM ubuntu:22.04
MAINTAINER Mingxun Wang "mwang87@gmail.com"

RUN apt-get update && apt-get install -y build-essential libarchive-dev wget vim git-core

# Install Mamba
ENV CONDA_DIR /opt/conda
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O ~/miniforge.sh && /bin/bash ~/miniforge.sh -b -p /opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH

# Adding to bashrc
RUN echo "export PATH=$CONDA_DIR:$PATH" >> ~/.bashrc

COPY requirements.txt environment-gemelli-worker.yml ./
RUN conda install -y python=3.12.2 && conda clean -afy
RUN pip install -r requirements.txt
RUN conda env create -f environment-gemelli-worker.yml && conda clean -afy
ENV GEMELLI_WORKER_PYTHON=/opt/conda/envs/gemelli-standalone/bin/python
RUN pip install git+https://github.com/Wang-Bioinformatics-Lab/GNPSDataPackage.git

COPY . /app
WORKDIR /app
