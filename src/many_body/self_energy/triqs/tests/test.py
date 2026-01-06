import src.many_body.self_energy.triqs.run as base
# Set configuration variables for test run
kmesh = [4, 4, 4]  # Example k-mesh values
mixing = 0.5    # Example mixing parameter
print(base.interaction)
base.interaction = "DMFT"


def test():
    # Main function call goes here
    print("Welcome to Testing! This is a Firefly run with k-mesh:", kmesh)
    result = base.run()
    return abs(result - 3.14) < 1e-6  # Example test condition

if __name__ == "__main__":
    test()
