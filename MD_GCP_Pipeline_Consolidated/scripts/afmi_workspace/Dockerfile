# AlphaFold-Metainference Docker Container
# Based on Ubuntu 20.04 with x86_64 architecture for OpenMM-Plumed-MPI compatibility

# Explicitly specify platform for x86_64 emulation on ARM Macs
FROM --platform=linux/amd64 continuumio/miniconda3:latest

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# Set working directory
WORKDIR /workspace

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    curl \
    vim \
    unzip \
    gfortran \
    && rm -rf /var/lib/apt/lists/*

# Conda is already installed in the base image
# Initialize conda for bash
RUN conda init bash && \
    echo "source /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc

# Create and activate the afmi-tutorial environment
RUN conda create -y -n afmi-tutorial python=3.9 && \
    echo "conda activate afmi-tutorial" >> ~/.bashrc

# Install mamba in base to speed up and reduce memory for solving
RUN /opt/conda/bin/conda install -y -n base -c conda-forge mamba

# Install conda packages in the environment
RUN /opt/conda/bin/mamba install -y -n afmi-tutorial -c conda-forge \
    openmpi \
    mpi4py \
    openmm=8.0 \
    plumed=2.8.2 \
    py-plumed \
    cmake \
    swig \
    pandas \
    mdtraj \
    biopython \
    matplotlib \
    gromacs \
    && /opt/conda/bin/mamba install -y -n afmi-tutorial -c anaconda ipykernel \
    && /opt/conda/bin/conda clean -afy

# Set up the shell to use bash with conda
SHELL ["/bin/bash", "-c"]

# Clone and build OpenMM-Plumed-MPI
RUN git clone https://github.com/vendruscolo-lab/OpenMM-Plumed-MPI /workspace/OpenMM-Plumed-MPI

WORKDIR /workspace/OpenMM-Plumed-MPI

# Extract and prepare OpenMM-Plumed-MPI
RUN mkdir -p build install openmm plumed/include plumed/lib && \
    unzip openmm.zip -d openmm && \
    unzip plumed_lib.zip -d plumed/lib && \
    unzip plumed_include.zip -d plumed/include

# Build OpenMM-Plumed-MPI
WORKDIR /workspace/OpenMM-Plumed-MPI/build

RUN source /opt/conda/etc/profile.d/conda.sh && \
    conda activate afmi-tutorial && \
    export CONDA_PREFIX=/opt/conda/envs/afmi-tutorial && \
    export PATH=$CONDA_PREFIX/bin:$PATH && \
    export CPLUS_INCLUDE_PATH=$CONDA_PREFIX/include:$CONDA_PREFIX/lib/python3.9/site-packages/mpi4py/include:$CPLUS_INCLUDE_PATH && \
    export C_INCLUDE_PATH=$CONDA_PREFIX/include:$CONDA_PREFIX/lib/python3.9/site-packages/mpi4py/include:$C_INCLUDE_PATH && \
    export LIBRARY_PATH=$CONDA_PREFIX/lib:$LIBRARY_PATH && \
    export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH && \
    export MPI_HOME=$CONDA_PREFIX && \
    export MPICC=$CONDA_PREFIX/bin/mpicc && \
    export MPICXX=$CONDA_PREFIX/bin/mpicxx && \
    cmake .. \
        -DCMAKE_PREFIX_PATH=$CONDA_PREFIX \
        -DMPI_C_COMPILER=$CONDA_PREFIX/bin/mpicc \
        -DMPI_CXX_COMPILER=$CONDA_PREFIX/bin/mpicxx \
        -DMPI_C_INCLUDE_PATH=$CONDA_PREFIX/include \
        -DMPI_CXX_INCLUDE_PATH=$CONDA_PREFIX/include && \
    make && \
    make install && \
    make PythonInstall && \
    cp -r ../install/lib/* /opt/conda/envs/afmi-tutorial/lib/

# Clone AlphaFold-IDP repository
WORKDIR /workspace
RUN git clone https://github.com/vendruscolo-lab/AlphaFold-IDP

# Copy the prepared files from host (will be mounted)
# This directory will be mounted when running the container

# Set the default working directory to prep_run
WORKDIR /workspace/tutorial

# Create a startup script
RUN echo '#!/bin/bash' > /workspace/start.sh && \
    echo 'source /opt/conda/etc/profile.d/conda.sh' >> /workspace/start.sh && \
    echo 'conda activate afmi-tutorial' >> /workspace/start.sh && \
    echo 'export CONDA_PREFIX=/opt/conda/envs/afmi-tutorial' >> /workspace/start.sh && \
    echo 'export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH' >> /workspace/start.sh && \
    echo 'export PYTHONPATH=$CONDA_PREFIX/lib/python3.9/site-packages:$PYTHONPATH' >> /workspace/start.sh && \
    echo 'cd /workspace/tutorial' >> /workspace/start.sh && \
    echo 'echo "=================================="' >> /workspace/start.sh && \
    echo 'echo "AlphaFold-Metainference Container"' >> /workspace/start.sh && \
    echo 'echo "=================================="' >> /workspace/start.sh && \
    echo 'echo ""' >> /workspace/start.sh && \
    echo 'echo "Environment: afmi-tutorial (activated)"' >> /workspace/start.sh && \
    echo 'echo "Working directory: /workspace/tutorial"' >> /workspace/start.sh && \
    echo 'echo ""' >> /workspace/start.sh && \
    echo 'echo "To run energy minimization:"' >> /workspace/start.sh && \
    echo 'echo "  python simulate_em.py 7.4 298"' >> /workspace/start.sh && \
    echo 'echo ""' >> /workspace/start.sh && \
    echo 'echo "To run AF-MI simulation (2 replicas):"' >> /workspace/start.sh && \
    echo 'echo "  mpirun -np 2 python simulate.py 7.4 298"' >> /workspace/start.sh && \
    echo 'echo ""' >> /workspace/start.sh && \
    echo 'echo "To run analysis:"' >> /workspace/start.sh && \
    echo 'echo "  cd analysis && ./script.sh 2"' >> /workspace/start.sh && \
    echo 'echo ""' >> /workspace/start.sh && \
    echo 'exec /bin/bash' >> /workspace/start.sh && \
    chmod +x /workspace/start.sh

# Expose port for potential Jupyter notebook (optional)
EXPOSE 8888

# Set environment variables
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
ENV CONDA_PREFIX=/opt/conda/envs/afmi-tutorial
ENV LD_LIBRARY_PATH=/opt/conda/envs/afmi-tutorial/lib:$LD_LIBRARY_PATH
ENV PYTHONPATH=/opt/conda/envs/afmi-tutorial/lib/python3.9/site-packages:$PYTHONPATH

# Default command
CMD ["/workspace/start.sh"]

