"""
Run the Python binding test suites and exit nonzero if any test fails.

Used by scripts/test-all.sh. Loads a defaults-only config first, so the
results don't depend on whatever build/bin/input.cfg currently holds.
"""
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "imports"))
sys.path.insert(0, HERE)

import cpp_imports as ff

cfg = tempfile.NamedTemporaryFile("w", suffix=".cfg", delete=False)
cfg.close()
try:
    ff.load_config(cfg.name)
finally:
    os.unlink(cfg.name)

import test_field_data_flow
import test_module_interface

failed = []

# test_field_data_flow has its own runner and summary
if not test_field_data_flow.run_all_tests():
    failed.append("test_field_data_flow")

# test_module_interface is a set of test_* functions returning True/False
names = sorted(n for n in dir(test_module_interface)
               if n.startswith("test_") and callable(getattr(test_module_interface, n)))
iface_failed = [n for n in names if not getattr(test_module_interface, n)()]
print(f"\ntest_module_interface: {len(names) - len(iface_failed)}/{len(names)} passed")
for n in iface_failed:
    print(f"  FAILED: {n}")
    failed.append(f"test_module_interface.{n}")

print()
if failed:
    print(f"Python tests FAILED: {len(failed)} failure(s)")
    sys.exit(1)
print("All Python tests passed")
