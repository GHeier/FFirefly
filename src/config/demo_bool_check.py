"""Demo module showing practical use of call_python_func_bool."""

def is_system_ready():
    """Check if system is ready for computation."""
    # In practice, this could check for files, resources, etc.
    return True

def has_converged():
    """Check if iterative calculation has converged."""
    # In practice, this would check convergence criteria
    tolerance = 1e-6
    current_error = 1e-7
    return current_error < tolerance

def should_continue_iteration():
    """Check if iteration should continue."""
    max_iterations = 100
    current_iteration = 50
    return current_iteration < max_iterations

def is_temperature_valid():
    """Check if temperature is in valid range."""
    temperature = 300  # Kelvin
    return 0 < temperature < 1000
