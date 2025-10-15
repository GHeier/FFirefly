"""Test module for testing bool call_python_func_bool interface."""

def return_true():
    """Returns True."""
    return True

def return_false():
    """Returns False."""
    return False

def return_truthy_int():
    """Returns a truthy integer (1)."""
    return 1

def return_falsy_int():
    """Returns a falsy integer (0)."""
    return 0

def return_string():
    """Returns a non-boolean value for testing error handling."""
    return "not a boolean"

def complex_check():
    """Performs a computation and returns bool."""
    x = 5 + 3
    return x > 7  # True

def another_complex_check():
    """Performs a computation and returns bool."""
    values = [1, 2, 3, 4, 5]
    return len(values) < 3  # False

if __name__ == "__main__":
    # Test the functions
    print(f"return_true(): {return_true()}")
    print(f"return_false(): {return_false()}")
    print(f"return_truthy_int(): {return_truthy_int()}")
    print(f"return_falsy_int(): {return_falsy_int()}")
    print(f"return_string(): {return_string()}")
    print(f"complex_check(): {complex_check()}")
    print(f"another_complex_check(): {another_complex_check()}")
