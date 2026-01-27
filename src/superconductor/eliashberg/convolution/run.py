import firefly
import firefly.config as cfg

import eliashberg
# Relevant variables loaded from the above import



def run():
    # Main function call goes here
    eig = eliashberg.run_lanczos()
    #eig = eliashberg.run_power_iteration()
    return eig  # Return something of any type that can be tested in the test suite.

if __name__ == "__main__": # Runs on file execution
    run()


