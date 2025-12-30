import firefly
import firefly.config as cfg

# Load relevant variables from the configuration
kmesh = cfg.k_mesh



def run():
    # Main function call goes here
    print("Hello, World! This is a Firefly run with k-mesh:", kmesh)
    return 3.14  # Return something of any type that can be tested in the test suite.

if __name__ == "__main__": # Runs on file execution
    run()


