#!/bin/bash
# Build FFirefly and run every test suite:
#   1. C++ unit tests          (fly.x --test)
#   2. Python binding tests    (src/module/tests/run_python_tests.py)
#   3. Julia binding tests     (src/module/tests/run_julia_tests.jl)
# Exits nonzero if any suite fails. Suites whose interpreter is not installed
# are reported as SKIPPED (or FAILED with --strict).
#
# Not yet included (stale against the current bindings API):
#   pypkg/firefly/tests/test_tensor_fields.py
#   jlpkg/Firefly/test/test_tensor_fields.jl
#   src/objects/CMField/tests/field_tests.jl

usage() {
    echo "Usage: $0 [--no-build] [--no-python] [--no-julia] [--strict]"
    echo "  --no-build   Use the existing build instead of rebuilding first"
    echo "  --no-python  Skip the Python binding tests"
    echo "  --no-julia   Skip the Julia binding tests"
    echo "  --strict     Treat a missing python3/julia as a failure"
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(realpath "${SCRIPT_DIR}/..")"

build=1
run_python=1
run_julia=1
strict=0
for arg in "$@"; do
    case "$arg" in
        --no-build) build=0 ;;
        --no-python) run_python=0 ;;
        --no-julia) run_julia=0 ;;
        --strict) strict=1 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $arg"; usage; exit 2 ;;
    esac
done

names=()
results=()
record() {
    names+=("$1")
    results+=("$2")
}

# Run a suite from the repo root with stdin closed, so nothing can block on it
run_suite() {
    local name="$1"
    shift
    echo
    echo "================================================================"
    echo "  $name"
    echo "================================================================"
    if (cd "$ROOT" && "$@" < /dev/null); then
        record "$name" PASSED
    else
        record "$name" FAILED
    fi
}

skip_suite() {
    if [ "$strict" -eq 1 ]; then
        record "$1" "FAILED ($2)"
    else
        record "$1" "SKIPPED ($2)"
    fi
}

if [ "$build" -eq 1 ]; then
    echo "Building..."
    if ! "$ROOT/scripts/fly-build.sh"; then
        echo "Build FAILED"
        exit 1
    fi
fi

if [ ! -x "$ROOT/build/bin/fly.x" ]; then
    echo "build/bin/fly.x not found - build first or drop --no-build"
    exit 1
fi

run_suite "C++ unit tests" "$ROOT/build/bin/fly.x" --test

if [ "$run_python" -eq 1 ]; then
    if command -v python3 > /dev/null; then
        run_suite "Python binding tests" python3 "$ROOT/src/module/tests/run_python_tests.py"
    else
        skip_suite "Python binding tests" "python3 not found"
    fi
fi

if [ "$run_julia" -eq 1 ]; then
    if command -v julia > /dev/null; then
        run_suite "Julia binding tests" julia --project="$ROOT/jlpkg/Firefly" \
            "$ROOT/src/module/tests/run_julia_tests.jl"
    else
        skip_suite "Julia binding tests" "julia not found"
    fi
fi

echo
echo "================================================================"
echo "  Summary"
echo "================================================================"
status=0
for i in "${!names[@]}"; do
    printf "  %-24s %s\n" "${names[$i]}" "${results[$i]}"
    case "${results[$i]}" in
        FAILED*) status=1 ;;
    esac
done
exit $status
