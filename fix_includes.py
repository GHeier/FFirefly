#!/usr/bin/env python3
"""
Fix relative includes in FFirefly codebase.

Changes relative includes like:
  #include "../objects/vec.hpp"
To absolute includes from project root:
  #include "src/objects/vec.hpp"
"""

import re
from pathlib import Path
import sys

def fix_relative_includes(file_path):
    """Fix relative includes in a single file."""
    with open(file_path, 'r') as f:
        content = f.read()

    original_content = content

    # Pattern to match relative includes like #include "../path/file.hpp"
    # Captures the full relative path
    pattern = r'#include\s+"(\.\./[^"]+)"'

    def replace_include(match):
        rel_path = match.group(1)

        # Count how many levels up we go
        parts = rel_path.split('/')
        up_count = sum(1 for p in parts if p == '..')

        # Remove the ../ parts to get the remaining path
        path_parts = [p for p in parts if p != '..']

        # Simple approach: ../ means we're going up from current file location
        # The number of ../ tells us how many directories to go up from the file's location
        # Then we append the remaining path

        # Example: file is src/hamiltonian/models/file.cpp
        #          include is "../config/load/c_config.h"
        #          up_count = 1, path_parts = ['config', 'load', 'c_config.h']
        #          We're in hamiltonian/models, go up 1 = hamiltonian, then config/load/c_config.h
        #          Result: src/config/load/c_config.h (wrong!)
        #          Actually: from hamiltonian/models, ../ goes to hamiltonian/, so ../config means hamiltonian/config (wrong)
        #          From hamiltonian/models, ../ goes to hamiltonian, so we need to add path from hamiltonian
        #          But ../config means config is a sibling of models, which is under hamiltonian
        #          So the correct interpretation: from src/hamiltonian/models, ../ = src/hamiltonian, ../config = src/config

        # Get file's directory relative to src/
        file_rel_to_src = file_path.relative_to(Path('src'))
        file_dir_parts = list(file_rel_to_src.parent.parts)

        # Apply the up traversal
        # If we're at ['hamiltonian', 'models'] and up_count=1, we get ['hamiltonian']
        # If up_count=2, we get []
        remaining_dirs = file_dir_parts[:-up_count] if up_count <= len(file_dir_parts) else []

        # Build final path
        final_parts = remaining_dirs + path_parts
        new_path = 'src/' + '/'.join(final_parts)

        return f'#include "{new_path}"'

    # Replace all matches
    new_content = re.sub(pattern, replace_include, content)

    if new_content != original_content:
        with open(file_path, 'w') as f:
            f.write(new_content)
        return True
    return False

def main():
    src_dir = Path('src')

    # Find all C/C++ source and header files
    files_to_fix = []
    for ext in ['*.cpp', '*.c', '*.hpp', '*.h']:
        for file_path in src_dir.rglob(ext):
            # Skip archive directories
            if '/archive/' in str(file_path) or '\\archive\\' in str(file_path):
                continue
            files_to_fix.append(file_path)

    print(f"Found {len(files_to_fix)} files to check...")

    fixed_count = 0
    for file_path in files_to_fix:
        try:
            if fix_relative_includes(file_path):
                print(f"  Fixed: {file_path}")
                fixed_count += 1
        except Exception as e:
            print(f"  Error processing {file_path}: {e}", file=sys.stderr)

    print(f"\n✓ Fixed {fixed_count} files")

if __name__ == '__main__':
    main()
