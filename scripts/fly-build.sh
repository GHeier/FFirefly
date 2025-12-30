#!/bin/bash

# Resolve the absolute path of the script, following symlinks if necessary
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
folder="${SCRIPT_DIR}/.."
folder="$(realpath "$folder")"

# Default to failure
status=1

if [ "$1" == "-v" ]; then
    cmake -S "$folder" -B "$folder/build" -G Ninja 1>/dev/null \
    && {
        output=$(cmake --build "$folder/build")
        build_status=$?
        echo "$output"
        exit $build_status
    } \
    && status=0

elif [ "$1" == "-vv" ]; then
    cmake -S "$folder" -B "$folder/build" -G Ninja \
    && {
        output=$(cmake --build "$folder/build")
        build_status=$?
        echo "$output"
        exit $build_status
    } \
    && status=0

elif [ -z "$1" ]; then
    cmake -S "$folder" -B "$folder/build" -G Ninja 1>/dev/null \
    && {
        output=$(cmake --build "$folder/build")
        build_status=$?
        echo "$output" | grep -v -E "Building|Generating|Linking"
        exit $build_status
    } \
    && status=0

else
    echo "Error: Invalid argument. Use -v for verbose or -vv for very verbose output."
    exit 1
fi

exit $status
