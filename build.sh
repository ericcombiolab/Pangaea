#!/bin/bash
set -e
set -o pipefail

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [PANGAEA-BUILD] $*"
}
log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [PANGAEA-BUILD][ERROR] $*" >&2
}

# >>> conda initialize >>>
if command -v conda >/dev/null 2>&1; then
    __conda_setup="$(conda 'shell.bash' 'hook' 2> /dev/null)"
    if [ $? -eq 0 ]; then
        eval "$__conda_setup"
    else
        if [ -f "$(dirname \"$(command -v conda)\")/../etc/profile.d/conda.sh" ]; then
            . "$(dirname \"$(command -v conda)\")/../etc/profile.d/conda.sh"
        fi
    fi
    unset __conda_setup
else
    log_error "Conda not found in PATH. Please install or load Conda."
    exit 1
fi
# <<< conda initialize <<<

# After conda initialization, the default environment is usually 'base'. Please ensure 'mamba' is installed in 'base', or activate another environment where 'mamba' is available before running this script.
# Example: conda activate <env_with_mamba>

# Create pangaea environment if not present
if ! conda env list | grep -q '^pangaea[[:space:]]'; then
    log "Creating conda environment 'pangaea'..."
    mamba env create -f environment.yaml -y
else
    log "Conda environment 'pangaea' already exists."
fi

# Create athena-meta environment if not present
if ! conda env list | grep -q '^athena-meta[[:space:]]'; then
    log "Creating conda environment 'athena-meta'..."
    mamba create -n athena-meta bioconda::athena_meta -y
else
    log "Conda environment 'athena-meta' already exists."
fi

log "Activating pangaea environment."
conda activate pangaea

# Install PyTorch (CPU version by default)
if ! python -c "import torch" 2>/dev/null; then
    log "Installing PyTorch (CPU version)..."
    pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
else
    log "PyTorch already installed."
fi

# Install rph_kmeans
log "Installing rph_kmeans..."
cd third_parties/rph_kmeans
python setup.py install
cd ../../

# Build C++ binaries
log "Building C++ tools..."
cd src/cpptools
if [ -d "build" ]; then
    log "Removing old build directory."
    rm -rf build
fi
mkdir build && cd build
cmake ..
make
cd ../../../

log "Pangaea is successfully built and installed under conda environment 'pangaea'."