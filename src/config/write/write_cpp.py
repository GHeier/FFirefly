from write_c import remove_lines_between_phrases, add_lines_between_phrases
start_phrase = '// Global Variables are listed below'
end_phrase = '// End of Global Variables'

start_func_phrase = '    // Load the C++ configuration file'
end_func_phrase = '    // End of Global Functions'

def format_var_line(key, value):
    if (type(value) == str):
        return f"string {key}(c_{key});"
    elif (type(value) == int):
        return f"int {key} = c_{key};"
    elif (type(value) == float):
        return f"float {key} = c_{key};"
    elif (type(value) == bool):
        return f"bool {key} = c_{key};"
    elif (type(value) == list):
        if (type(value[0]) == int):
            return f"int {key}[3] = c_{key};"
        elif (type(value[0]) == float):
            return f"float {key}[3] = c_{key};"
        elif (type(value[0]) == list):
            if (type(value[0][0]) == int):
                return f"int {key}[3][3] = c_{key};"
            elif (type(value[0][0]) == float):
                return f"float {key}[3][3] = c_{key};"
        else:
            print("Error: Unsupported type in config file (list section)")
            exit(1)

def format_header_line(key, value):
    if (type(value) == str):
        return f"extern string {key};"
    elif (type(value) == int):
        return f"extern int {key};"
    elif (type(value) == float):
        return f"extern float {key};"
    elif (type(value) == bool):
        return f"extern bool {key};"
    elif (type(value) == list):
        if (type(value[0]) == int):
            return f"extern int {key}[3];"
        elif (type(value[0]) == float):
            return f"extern float {key}[3];"
        elif (type(value[0]) == list):
            if (type(value[0][0]) == int):
                return f"extern int {key}[3][3];"
            elif (type(value[0][0]) == float):
                return f"extern float {key}[3][3];"
        else:
            print("Error: Unsupported type in config file (list section)")
            exit(1)

def format_func_line(key, value):
    if (type(value) == list):
        return f"    for (int i = 0; i < 3; i++) {key}[i] = c_{key}[i];"
    return f"    {key} = c_{key};"

def get_new_func_lines(ALL):
    new_lines = []
    for section in ALL:
        new_lines.append('\n    //' + '[' + section + ']\n')
        for key, value in ALL[section].items():
            if (section == "CELL"): 
                continue
            new_lines.append(format_func_line(key, value) + '\n')
        new_lines.append('')
    return new_lines
        

def get_new_var_lines(ALL):
    new_lines = []
    for section in ALL:
        new_lines.append('\n//' + '[' + section + ']\n')
        for key, value in ALL[section].items():
            new_lines.append(format_header_line(key, value) + '\n')
        new_lines.append('')
    return new_lines

def get_new_header_lines(ALL):
    new_lines = []
    for section in ALL:
        new_lines.append('\n//' + '[' + section + ']\n')
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

def add_lines_to_header(ALL, file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    lines = remove_lines_between_phrases(lines, start_phrase, end_phrase)
    new_lines = get_new_var_lines(ALL)
    updated_lines = add_lines_between_phrases(lines, new_lines, start_phrase, end_phrase)
    
    # Write the modified content back to the file
    with open(file_path, 'w') as file:
        file.writelines(updated_lines)

def write_cpp_header(ALL):
    file_path = 'load/cpp_config.h'
    add_lines_to_header(ALL, file_path)
    print(f"Successfully updated the file '{file_path}'.")

def write_cpp(ALL):
    file_path = 'load/cpp_config.cpp'
    add_lines_to_vars(ALL, file_path)
    add_lines_to_funcs(ALL, file_path)
    print(f"Successfully updated the file '{file_path}'.")
