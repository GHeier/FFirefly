#!/bin/bash

echo "========================================"
echo "Integration Test: method2 functions"
echo "========================================"

cd /home/g/Research/FFirefly

# Test 1: Python method2
echo ""
echo "--- Test 1: run_python_method2 ---"
echo "Running: python3 src/hamiltonian/DOS/gaussian/run.py < test_method2/test_gaussian.cfg"
python3 src/hamiltonian/DOS/gaussian/run.py < test_method2/test_gaussian.cfg
PYTHON_EXIT=$?
echo "Exit code: $PYTHON_EXIT"

if [ $PYTHON_EXIT -eq 0 ]; then
    echo "✓ Python method2 test PASSED"
else
    echo "✗ Python method2 test FAILED"
fi

# Test 2: Julia method2 (simple test)
echo ""
echo "--- Test 2: run_julia_method2 (simple) ---"
echo "Running: julia test_method2/test_julia_simple.jl < test_method2/test_gaussian.cfg"
julia test_method2/test_julia_simple.jl < test_method2/test_gaussian.cfg
JULIA_EXIT=$?
echo "Exit code: $JULIA_EXIT"

if [ $JULIA_EXIT -eq 0 ]; then
    echo "✓ Julia method2 test PASSED"
else
    echo "✗ Julia method2 test FAILED"
fi

# Summary
echo ""
echo "========================================"
echo "Summary:"
echo "========================================"
if [ $PYTHON_EXIT -eq 0 ]; then
    echo "  Python method2: PASSED ✓"
else
    echo "  Python method2: FAILED ✗"
fi

if [ $JULIA_EXIT -eq 0 ]; then
    echo "  Julia method2:  PASSED ✓"
else
    echo "  Julia method2:  FAILED ✗"
fi
echo "========================================"

# Exit with success if both passed
if [ $PYTHON_EXIT -eq 0 ] && [ $JULIA_EXIT -eq 0 ]; then
    exit 0
else
    exit 1
fi
