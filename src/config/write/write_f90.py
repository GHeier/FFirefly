from .write_c import add_lines
start_phrase = '    ! Global variables'
end_phrase = '    ! End of global variables'
start_func_phrase = '    ! Global functions'
end_func_phrase = '    ! End of global functions'
start_load_phrase = '        ! Load variables'
end_load_phrase = '        ! End of loading variables'


def format_var_line(key, value, section):
    index = ''
    if section == 'BANDS':
        index = '(50)'
        if key == 'band':
            return f"    character(len=50) :: band(50,50)"
    if (type(value) == str):
        return f"    character(len=50) :: {key}{index}"
    elif (type(value) == int):
        return f"    integer(c_int), bind(C, name=\"c_{key}\") :: c_{key}{index}\n    integer :: {key}{index}"
    elif (type(value) == float):
        return f"    real(c_float), bind(C, name=\"c_{key}\") :: c_{key}{index}\n    real :: {key}{index}"
    elif (type(value) == bool):
        return f"    logical(c_bool), bind(C, name=\"c_{key}\") :: c_{key}{index}\n    logical :: {key}{index}"
    elif (type(value) == list):
        if (type(value[0]) == int):
            return f"    integer(c_int), bind(C, name=\"c_{key}\") :: c_{key}{index}(3)\n    integer :: {key}{index}(3)"
        elif (type(value[0]) == float):
            return f"    real(c_float), bind(C, name=\"c_{key}\") :: c_{key}{index}(3)\n    real :: {key}{index}(3)"
        elif (type(value[0]) == list):
            if (type(value[0][0]) == int):
                return f"    integer(c_int), bind(C, name=\"c_{key}\") :: c_{key}{index}(3,3)\n    integer :: {key}{index}(3,3)"
            elif (type(value[0][0]) == float):
                return f"    real(c_float), bind(C, name=\"c_{key}\") :: c_{key}{index}(3,3)\n    real :: {key}{index}(3,3)"
        else:
            print("Error: Unsupported type in config file (list section)")
            exit(1)
    else:
        print("Error: Unsupported type in config file")
        exit(1)

def format_func_line(key, value, section):
    if (type(value) == str):
        return f"        function get_{key}() bind(C)\n            use iso_c_binding\n            type(c_ptr) :: get_{key}\n    end function get_{key}"
    return ""

def format_load_line(key, value, section):
    if (type(value) == str):
        return f"        {key} = get_string(get_{key}())"
    else:
        return f"        {key} = c_{key}"

def replace_double_slash(file_path):
    # Read the file content
    with open(file_path, 'r') as file:
        content = file.read()
    
    # Replace all occurrences of "//" with "!"
    modified_content = content.replace("//", "!")
    
    
    # Write the modified content to the output file
    with open(file_path, 'w') as file:
        file.write(modified_content)
    

def write_f90(ALL):
    file_path = 'load/fortran_config.f90'
    add_lines(ALL, file_path, start_phrase, end_phrase, format_var_line)
    add_lines(ALL, file_path, start_func_phrase, end_func_phrase, format_func_line)
    add_lines(ALL, file_path, start_load_phrase, end_load_phrase, format_load_line)
    replace_double_slash(file_path)

    print(f"Successfully updated the file '{file_path}'.")

