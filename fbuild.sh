#!/bin/bash

# Determine the project folder (default to current directory if not set)
folder=/home/g/Research/bcs_diagonalization

# Check if the first argument is provided and handle verbosity levels
if [ "$1" == "-v" ]; then
    cmake -S "$folder" -B "$folder/build" > /dev/null >&1 && cmake --build "$folder/build" --parallel 4
elif [ "$1" == "-vv" ]; then
    cmake -S "$folder" -B "$folder/build" && cmake --build "$folder/build" --parallel 4
elif [ -z "$1" ]; then
    cmake -S "$folder" -B "$folder/build" > /dev/null >&1 && cmake --build "$folder/build" --parallel 4 > /dev/null >&1
else
    echo "Error: Invalid argument. Use -v for verbose or -vv for very verbose output."
    exit 1
fi
