#!/bin/bash

set -e

OS="$(uname -s)"

if [ "$OS" = "Linux" ]; then
    echo "Detected Linux — using apt..."
    sudo apt update
    sudo apt install -y gcc g++ cmake libopenblas-dev liblapack-dev liblapacke-dev ninja-build libomp-dev ccache libboost-all-dev pybind11-dev libhdf5-dev

elif [ "$OS" = "Darwin" ]; then
    echo "Detected macOS — using Homebrew..."

    if ! command -v brew &>/dev/null; then
        echo "Homebrew not found. Installing Homebrew first..."
        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    fi

    /opt/homebrew/bin/brew install gcc cmake openblas lapack ninja libomp ccache boost pybind11 hdf5

else
    echo "Unsupported OS: $OS"
    exit 1
fi

echo "All packages installed successfully!"
