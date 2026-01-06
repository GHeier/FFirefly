import src.superconductor.eliashberg.lanczos.run as base
# Set configuration variables for test run
base.kmesh = [4, 4, 4]  # Example k-mesh values



def test():
    # Main function call goes here
    print("Welcome to Testing! This is a Firefly run with k-mesh:", base.kmesh)
    result = base.run()
    return abs(result - 3.14) < 1e-6  # Example test condition

