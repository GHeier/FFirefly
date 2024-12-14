from .write_c import add_lines
def format_var_line(key, value, section):
    if (key == "band"):
        return f"band = [''] * 50\n"

    if (type(value) == str):
        if (key == "band"):
            return f"{key} = [''] * 50"
        return f"{key} = ''"
    elif (type(value) == int):
        if (key == "band"):
            return f"{key} = [0] * 50"
        return f"{key} = 0"
    elif (type(value) == float):
        if (key == "band"):
            return f"{key} = [0.0] * 50"
        return f"{key} = 0.0"
    elif (type(value) == bool):
        if (key == "band"):
            return f"{key} = [False] * 50"
        return f"{key} = False"
    elif (type(value) == list):
        if (type(value[0]) == int):
            if (key == "band"):
                return f"{key} = [[0] * 3] * 50"
            return f"{key} = [0] * 3"
        elif (type(value[0]) == float):
            if (key == "band"):
                return f"{key} = [[0.0] * 3] * 50"
            return f"{key} = [0.0] * 3"
        elif (type(value[0]) == list):
            if (type(value[0][0]) == int):
                if (key == "band"):
                    return f"{key} = [[[0] * 3] * 3] * 50"
                return f"{key} = [[0] * 3] * 3"
            elif (type(value[0][0]) == float):
                if (key == "band"):
                    return f"{key} = [[[0.0] * 3] * 3] * 50"
                return f"{key} = [[0.0] * 3] * 3"
        else:
            print("Error: Unsupported type in config file ")
            exit(1)

def format_func_line(key, value, section):
    base_val = f"global {key}\n"
    types_val = ""
    if type(value) == str:
        if section == 'BANDS':
            types_val = f"{key} = ctypes.c_char_p.in_dll(lib, 'c_{key}').value.decode('utf-8')"
        types_val = f"{key} = ctypes.c_char_p.in_dll(lib, 'c_{key}').value.decode('utf-8')"
                
    elif type(value) == int:
        if section == 'BANDS':
            types_val = f"{key} = list((ctypes.c_int * 50).in_dll(lib, 'c_{key}'))"
        types_val = f"{key} = ctypes.c_int.in_dll(lib, 'c_{key}').value"
    elif type(value) == float:
        if section == 'BANDS':
            types_val = f"{key} = list((ctypes.c_float * 50).in_dll(lib, 'c_{key}'))"
        types_val = f"{key} = ctypes.c_float.in_dll(lib, 'c_{key}').value"
    elif type(value) == bool:
        if section == 'BANDS':
            types_val = f"{key} = list((ctypes.c_bool * 50).in_dll(lib, 'c_{key}'))"
        types_val = f"{key} = ctypes.c_bool.in_dll(lib, 'c_{key}').value"
    elif type(value) == list:
        if type(value[0]) == int:
            if section == 'BANDS':
                types_val = f"{key} = list((ctypes.c_int * 50).in_dll(lib, 'c_{key}'))"
            types_val = f"{key} = list((ctypes.c_int * 3).in_dll(lib, 'c_{key}'))"
        elif type(value[0]) == float:
            if section == 'BANDS':
                types_val = f"{key} = list((ctypes.c_float * 50).in_dll(lib, 'c_{key}'))"
            types_val = f"{key} = list((ctypes.c_float * 3).in_dll(lib, 'c_{key}'))"
        elif type(value[0]) == list:
            if type(value[0][0]) == int:
                if section == 'BANDS':
                    types_val = f"{key} = [[[((((ctypes.c_int * 3) * 3) * 50).in_dll(lib, 'c_{key}'))[i][j] for j in range(3)] for i in range(3)] for k in range(50)]"
                types_val = f"{key} = [[(((ctypes.c_int * 3) * 3).in_dll(lib, 'c_{key}'))[i][j] for j in range(3)] for i in range(3)]"
            elif type(value[0][0]) == float:
                if section == 'BANDS':
                    types_val = f"{key} = [[[((((ctypes.c_float * 3) * 3) * 50).in_dll(lib, 'c_{key}'))[i][j] for j in range(3)] for i in range(3)] for k in range(50)]"
                types_val = f"{key} = [[(((ctypes.c_float * 3) * 3).in_dll(lib, 'c_{key}'))[i][j] for j in range(3)] for i in range(3)]"
        else:
            print("Error: Unsupported type in config file ")
            exit(1)
    return base_val + types_val

def format_func2(key, value, section):
    return f"    m.attr(\"{key}\") = &{key};"

def replace_comments(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    with open(file_path, 'w') as file:
        for line in lines:
            file.write(replace_comment(line))

def replace_comment(line):
    if line.startswith('//'):
        return line.replace('//', '#', 1)
    return line

start_phrase = '# Begin Global Variables'
end_phrase = '# End Global Variables'
start_func_phrase = '# Begin Function Declarations'
end_func_phrase = '# End Function Declarations'
start_func2_phrase = '    // Begin PyBind Definitions'
end_func2_phrase = '    // End PyBind Definitions'

def write_py(ALL):
    file_path = 'load/python_config.py'
    #add_lines(ALL, file_path, start_phrase, end_phrase, format_var_line)
    add_lines(ALL, file_path, start_func_phrase, end_func_phrase, format_func_line)
    replace_comments(file_path)
    print(f"Successfully updated the file '{file_path}'.")
