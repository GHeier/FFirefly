#!/bin/bash

set -e

OS="$(uname -s)"

if [ "$OS" = "Linux" ]; then
    echo "Detected Linux — using apt..."
    sudo apt update
    sudo apt install -y gcc g++ cmake libopenblas-dev liblapack-dev liblapacke-dev ninja-build libomp-dev ccache libboost-all-dev pybind11-dev libhdf5-dev libfftw3-dev libopenblas-dev libcairo2-dev
    
    # Make Debian back-compatible with AUR install
    sudo mkdir -p /usr/include/openblas && cd /usr/include/openblas && sudo ln -sf ../lapacke.h .
    cd /usr/include && sudo ln -sf hdf5/serial/*.h .

    [ -d H2Lib ] || sudo git clone https://github.com/H2Lib/H2Lib.git && sudo cp H2Lib/Library/*.h /usr/local/include/
    cd H2Lib
    sudo make clean
    sudo make CFLAGS="-fPIC"

    sudo ln -sf Library/hmatrix.h .

    #TMPDIR=$(mktemp -d)
    #git clone https://github.com/H2Lib/H2Lib.git "$TMPDIR/H2Lib"
    #sudo cp -r "$TMPDIR/H2Lib/Library/"*.h /usr/local/include/
    ## build the library itself if needed
    #cd "$TMPDIR/H2Lib" && make
    #sudo cp libH2.a /usr/local/lib/    # or whatever the output artifact is called
    #rm -rf "$TMPDIR"

    # Setup Hierarchical Matrices
    #[ -d H2Lib ] || git clone https://github.com/H2Lib/H2Lib.git /home/g/FFirefly/external/H2Lib

elif [ "$OS" = "Darwin" ]; then
    echo "Detected macOS — using Homebrew..."

    BREW=/opt/homebrew/bin/brew

    if ! command -v "$BREW" &>/dev/null; then
        echo "Homebrew not found. Installing Homebrew first..."
        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    fi

    $BREW install gcc cmake openblas lapack ninja libomp ccache boost pybind11 hdf5

else
    echo "Unsupported OS: $OS"
    exit 1
fi

echo "All packages installed successfully!"
