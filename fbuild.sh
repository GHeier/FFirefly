#!/bin/bash

# Determine the project folder (default to current directory if not set)
folder=/home/g/Research/bcs_diagonalization

# Check if the first argument is provided and handle verbosity levels
if [ "$1" == "-v" ]; then
    echo "Running CMake with verbose output (level 1)..."
    cmake -S "$folder" -B "$folder/build" > /dev/null 2>&1 && cmake --build "$folder/build" --parallel 4
elif [ "$1" == "-vv" ]; then
    echo "Running CMake with very verbose output (level 2)..."
    cmake -S "$folder" -B "$folder/build" && cmake --build "$folder/build" --parallel 4
elif [ -z "$1" ]; then
    echo "Running CMake with minimal output (default)..."
    cmake -S "$folder" -B "$folder/build" > /dev/null 2>&1 && cmake --build "$folder/build" --parallel 4 > /dev/null 2>&1
else
    echo "Error: Invalid argument. Use -v for verbose or -vv for very verbose output."
    exit 1
fi
