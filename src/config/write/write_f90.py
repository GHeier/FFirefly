from write_c import remove_lines_between_phrases, add_lines_between_phrases
start_phrase = '    ! Global variables'
end_phrase = '    ! End of global variables'
start_func_phrase = '    ! Global functions'
end_func_phrase = '    ! End of global functions'
start_load_phrase = '        ! Load variables'
end_load_phrase = '        ! End of loading variables'


def format_var_line(key, value):
    if (type(value) == str):
        return f"    character(kind=c_char), bind(C, name=\"c_{key}\") :: c_{key}(50)\n    character(len=50) :: {key}"
    elif (type(value) == int):
        return f"    integer(c_int), bind(C, name=\"c_{key}\") :: c_{key}\n    integer :: {key}"
    elif (type(value) == float):
        return f"    real(c_float), bind(C, name=\"c_{key}\") :: c_{key}\n    real :: {key}"
    elif (type(value) == bool):
        return f"    logical(c_bool), bind(C, name=\"c_{key}\") :: c_{key}\n    logical :: {key}"
    elif (type(value) == list):
        if (type(value[0]) == int):
            return f"    integer(c_int), bind(C, name=\"c_{key}\") :: c_{key}(3)\n    integer :: {key}(3)"
        elif (type(value[0]) == float):
            return f"    real(c_float), bind(C, name=\"c_{key}\") :: c_{key}(3)\n    real :: {key}(3)"
        elif (type(value[0]) == list):
            if (type(value[0][0]) == int):
                return f"    integer(c_int), bind(C, name=\"c_{key}\") :: c_{key}(3,3)\n    integer :: {key}(3,3)"
            elif (type(value[0][0]) == float):
                return f"    real(c_float), bind(C, name=\"c_{key}\") :: c_{key}(3,3)\n    real :: {key}(3,3)"
        else:
            print("Error: Unsupported type in config file (list section)")
            exit(1)
    else:
        print("Error: Unsupported type in config file")
        exit(1)

def format_func_line(key, value):
    return f"        function get_{key}() bind(C)\n            use iso_c_binding\n            type(c_ptr) :: get_{key}\n        end function get_{key}\n"

def format_load_line(key, value):
    if (type(value) == str):
        return f"        {key} = get_string(get_{key}())"
    else:
        return f"        {key} = c_{key}"

def get_new_func_lines(ALL):
    new_lines = []
    for section in ALL:
        new_lines.append('\n    !' + '[' + section + ']\n')
        for key, value in ALL[section].items():
            if (section == "CELL" or type(value) != str): 
                continue
            new_lines.append(format_func_line(key, value) + '\n')
        new_lines.append('')
    return new_lines

def get_new_load_lines(ALL):
    new_lines = []
    for section in ALL:
        new_lines.append('\n        !' + '[' + section + ']\n')
        for key, value in ALL[section].items():
            new_lines.append(format_load_line(key, value) + '\n')
        new_lines.append('')
    return new_lines

def get_new_var_lines(ALL):
    new_lines = []
    for section in ALL:
        new_lines.append('\n    !' + '[' + section + ']\n')
        for key, value in ALL[section].items():
            new_lines.append(format_var_line(key, value) + '\n')
        new_lines.append('')
    return new_lines

def add_lines_to_vars(ALL, file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    lines = remove_lines_between_phrases(lines, start_phrase, end_phrase)
    new_lines = get_new_var_lines(ALL)
    updated_lines = add_lines_between_phrases(lines, new_lines, start_phrase, end_phrase)
    
    # Write the modified content back to the file
    with open(file_path, 'w') as file:
        file.writelines(updated_lines)

def add_lines_to_funcs(ALL, file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    lines = remove_lines_between_phrases(lines, start_func_phrase, end_func_phrase)
    new_lines = get_new_func_lines(ALL)
    updated_lines = add_lines_between_phrases(lines, new_lines, start_func_phrase, end_func_phrase)
    
    # Write the modified content back to the file
    with open(file_path, 'w') as file:
        file.writelines(updated_lines)

def add_lines_to_load(ALL, file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    lines = remove_lines_between_phrases(lines, start_load_phrase, end_load_phrase)
    new_lines = get_new_load_lines(ALL)
    updated_lines = add_lines_between_phrases(lines, new_lines, start_load_phrase, end_load_phrase)
    
    # Write the modified content back to the file
    with open(file_path, 'w') as file:
        file.writelines(updated_lines)

def write_f90(ALL):
    file_path = 'load/fortran_config.f90'
    add_lines_to_vars(ALL, file_path)
    add_lines_to_funcs(ALL, file_path)
    add_lines_to_load(ALL, file_path)

    print(f"Successfully updated the file '{file_path}'.")

