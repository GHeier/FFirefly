#!/usr/bin/env python3
"""
Simple test script for run_python_method2.
Reads from stdin and prints confirmation.
"""
import sys

print("=== Python method2 test script ===")
print("Reading from stdin...")

# Read all input
config_data = sys.stdin.read()

if config_data:
    print(f"Received {len(config_data)} bytes of config data")
    print("First 100 chars:")
    print(config_data[:100])
else:
    print("No input received")

print("Python script completed successfully")
sys.exit(0)
