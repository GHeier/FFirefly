from write_c import add_lines
def format_var_line(key, value, section):
    index = ''
    if section == 'BANDS':
        index = '[50]'
        if (key == "band"):
            return (
                    f"string {key}[50];\n"
                    )
    if (type(value) == str):
        return f"string {key}{index};"
    elif (type(value) == int):
        return f"int {key}{index};"
    elif (type(value) == float):
        return f"float {key}{index};"
    elif (type(value) == bool):
        return f"bool {key}{index};"
    elif (type(value) == list):
        if (type(value[0]) == int):
            return f"int {key}{index}[3];"
        elif (type(value[0]) == float):
            return f"float {key}{index}[3];"
        elif (type(value[0]) == list):
            if (type(value[0][0]) == int):
                return f"int {key}{index}[3][3];"
            elif (type(value[0][0]) == float):
                return f"float {key}{index}[3][3];"
        else:
            print("Error: Unsupported type in config file ")
            exit(1)

def format_header_line(key, value, section):
    index = ''
    if section == 'BANDS':
        index = '[50]'
        if (key == "band"):
            return (
                    f"extern string {key}[50];\n"
                    )
    if (type(value) == str):
        return f"extern string {key}{index};"
    elif (type(value) == int):
        return f"extern int {key}{index};"
    elif (type(value) == float):
        return f"extern float {key}{index};"
    elif (type(value) == bool):
        return f"extern bool {key}{index};"
    elif (type(value) == list):
        if (type(value[0]) == int):
            return f"extern int {key}{index}[3];"
        elif (type(value[0]) == float):
            return f"extern float {key}{index}[3];"
        elif (type(value[0]) == list):
            if (type(value[0][0]) == int):
                return f"extern int {key}{index}[3][3];"
            elif (type(value[0][0]) == float):
                return f"extern float {key}{index}[3][3];"
        else:
            print("Error: Unsupported type in config file (list section)")
            exit(1)

def format_func_line(key, value, section):
    if (type(value) == list):
        if (type(value[0]) == list):
            return f"    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {key}[i][j] = c_{key}[i][j];"
        return f"    for (int i = 0; i < 3; i++) {key}[i] = c_{key}[i];"
    elif (section == 'BANDS'):
        return f"    for (int i = 0; i < nbnd; i++) {key}[i] = c_{key}[i];"
    return f"    {key} = c_{key};"

start_phrase = '// Global Variables are listed below'
end_phrase = '// End of Global Variables'

start_func_phrase = '    // Load the C++ configuration file'
end_func_phrase = '    // End of Global Functions'


def write_cpp_header(ALL):
    file_path = 'load/cpp_config.h'
    add_lines(ALL, file_path, start_phrase, end_phrase, format_header_line)
    print(f"Successfully updated the file '{file_path}'.")

def write_cpp(ALL):
    file_path = 'load/cpp_config.cpp'
    add_lines(ALL, file_path, start_phrase, end_phrase, format_var_line)
    add_lines(ALL, file_path, start_func_phrase, end_func_phrase, format_func_line)
    print(f"Successfully updated the file '{file_path}'.")
