"""Test module for typed interface functions."""

# Integer tests
def return_int_42():
    """Returns integer 42."""
    return 42

def return_int_negative():
    """Returns negative integer."""
    return -100

def return_int_zero():
    """Returns zero."""
    return 0

# Float tests
def return_float_pi():
    """Returns pi as float."""
    return 3.14159

def return_float_negative():
    """Returns negative float."""
    return -2.71828

def return_float_zero():
    """Returns 0.0."""
    return 0.0

# Double tests (Python doesn't distinguish, but C does)
def return_double_large():
    """Returns large double value."""
    return 1.23456789012345

def return_double_small():
    """Returns small double value."""
    return 0.00000123456

# String tests
def return_string_hello():
    """Returns 'Hello, World!'."""
    return "Hello, World!"

def return_string_empty():
    """Returns empty string."""
    return ""

def return_string_special():
    """Returns string with special characters."""
    return "Test!@#$%^&*()"

# Computed values
def compute_sum():
    """Computes and returns sum."""
    return 10 + 20 + 30

def compute_product():
    """Computes and returns product as float."""
    return 2.5 * 4.0

def compute_message():
    """Constructs and returns a message."""
    name = "FFirefly"
    version = "1.0"
    return f"{name} v{version}"

if __name__ == "__main__":
    # Test all functions
    print(f"return_int_42(): {return_int_42()}")
    print(f"return_int_negative(): {return_int_negative()}")
    print(f"return_int_zero(): {return_int_zero()}")
    print(f"return_float_pi(): {return_float_pi()}")
    print(f"return_float_negative(): {return_float_negative()}")
    print(f"return_float_zero(): {return_float_zero()}")
    print(f"return_double_large(): {return_double_large()}")
    print(f"return_double_small(): {return_double_small()}")
    print(f"return_string_hello(): {return_string_hello()}")
    print(f"return_string_empty(): '{return_string_empty()}'")
    print(f"return_string_special(): {return_string_special()}")
    print(f"compute_sum(): {compute_sum()}")
    print(f"compute_product(): {compute_product()}")
    print(f"compute_message(): {compute_message()}")
