import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
import run as base
#import src.many_body.renormalization.FS_approx.run as base

# Set configuration variables for test run
base.kmesh = [4, 4, 4]  # Example k-mesh values



def test():
    # Main function call goes here
    print("Welcome to Testing! This is a Firefly run with k-mesh:", base.kmesh)
    expected = 3.14
    result = base.run()
    pass = abs(result - expected) < 1e-6  # Example test condition

    print("Expected: ", expected)
    print("Result: ", result)

    return pass

if __name__ == "__main__":
    pass = test()
    if pass:
        print("Test passed!")
    else:
        print("Test failed.")

