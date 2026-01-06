import firefly
import firefly.config as cfg

import src.superconductor.shared.bcs as bcs
# Relevant variables loaded in the above import

def run():
    # Main function call goes here
    eig = run_lanczos()
    return eig  # Return something of any type that can be tested in the test suite.

if __name__ == "__main__": # Runs on file execution
    run()


