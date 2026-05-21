#import src.many_body.self_energy.triqs.run as base
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
import run as base
# Set configuration variables for test run
kmesh = [4, 4, 4]  # Example k-mesh values
#base.mixing = 1.0    # Example mixing parameter
base.interaction = "DMFT"
base.mu = 0.0
base.n = 1.0
#base.U = 2.0
base.beta = 20.0
base.Temperature = 0.05

def test():
    # Main function call goes here
    print("Welcome to Testing! This is a Firefly run with k-mesh:", kmesh)
    result = base.run()
    return abs(result - 0.7397) < 1e-4  # Example test condition

if __name__ == "__main__":
    result = test()
    print(int(result) * "Test passed! Bethe Lattice condition met")
    print(int(not result) * "Test failed. Bethe Lattice condition unmet")
