#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
folder="${SCRIPT_DIR}/.."
folder="$(realpath "$folder")"

# Check if the first argument is provided and handle verbosity levels
if [ "$1" == "-v" ]; then
    cmake -S "$folder" -B "$folder/build" -G Ninja 1>/dev/null && cmake --build "$folder/build"
elif [ "$1" == "-vv" ]; then
    cmake -S "$folder" -B "$folder/build" -G Ninja && cmake --build "$folder/build"
elif [ -z "$1" ]; then
    cmake -S "$folder" -B "$folder/build" -G Ninja 1>/dev/null && cmake --build "$folder/build" | grep -v -E "Building|Generating|Linking"
else
    echo "Error: Invalid argument. Use -v for verbose or -vv for very verbose output."
    exit 1
fi
