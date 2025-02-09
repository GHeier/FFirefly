from .write_c import add_lines
from .write_py import replace_comments
def format_var_line(key, value, section):
    makevec = False
    if section == 'BANDS':
        makevec = True
        if (key == "band"):
            return (
                    f"{key}::Vector{{String}} = cfg.{key}\n"
                    )
    if (type(value) == str):
        if makevec:
            return f"{key}::Vector{{String}} = cfg.{key}"
        return f"{key}::String = cfg.{key}"
    elif (type(value) == int):
        if makevec:
            return f"{key}::Vector{{Int}} = cfg.{key}"
        return f"{key}::Int = cfg.{key}"
    elif (type(value) == float):
        if makevec:
            return f"{key}::Vector{{Float64}} = cfg.{key}"
        return f"{key}::Float64 = cfg.{key}"
    elif (type(value) == bool):
        if makevec:
            return f"{key}::Vector{{Bool}} = cfg.{key}"
        return f"{key}::Bool = cfg.{key}"
    elif (type(value) == list):
        if (type(value[0]) == int):
            return f"{key}::Array{{Int}} = cfg.{key}"
        elif (type(value[0]) == float):
            return f"{key}::Array{{Float64}} = cfg.{key}"
        elif (type(value[0]) == list):
            if (type(value[0][0]) == int):
                return f"{key}::Array{{Int}} = cfg.{key}"
            elif (type(value[0][0]) == float):
                return f"{key}::Array{{Float64}} = cfg.{key}"
        else:
            print("Error: Unsupported type in config file ")
            exit(1)

def format_config_line(key, value, section):
    makevec = False
    if section == 'BANDS':
        makevec = True
        if (key == "band"):
            return (
                    f"    vector<string> {key};\n"
                    )
    if (type(value) == str):
        if makevec:
            return f"    vector<string> {key};"
        return f"    string {key};"
    elif (type(value) == int):
        if makevec:
            return f"    vector<int> {key};"
        return f"    int {key};"
    elif (type(value) == float):
        if makevec:
            return f"    vector<float> {key};"
        return f"    float {key};"
    elif (type(value) == bool):
        if makevec:
            return f"    vector<bool> {key};"
        return f"    bool {key};"
    elif (type(value) == list):
        if (type(value[0]) == int):
            return f"    vector<int> {key};"
        elif (type(value[0]) == float):
            return f"    vector<float> {key};"
        elif (type(value[0]) == list):
            if (type(value[0][0]) == int):
                return f"    vector<vector<int>> {key};"
            elif (type(value[0][0]) == float):
                return f"    vector<vector<float>> {key};"
        else:
            print("Error: Unsupported type in config file ")
            exit(1)

def format_header_line(key, value, section):
    makevec = False
    if section == 'BANDS':
        makevec = True
        if (key == "band"):
            return (
                    f"extern vector<string> {key};\n"
                    )
    if (type(value) == str):
        if makevec:
            return f"extern vector<string> {key};"
        return f"extern string {key};"
    elif (type(value) == int):
        if makevec:
            return f"extern vector<int> {key};"
        return f"extern int {key};"
    elif (type(value) == float):
        if makevec:
            return f"extern vector<float> {key};"
        return f"extern float {key};"
    elif (type(value) == bool):
        if makevec:
            return f"extern vector<bool> {key};"
        return f"extern bool {key};"
    elif (type(value) == list):
        if (type(value[0]) == int):
            return f"extern vector<int> {key};"
        elif (type(value[0]) == float):
            return f"extern vector<float> {key};"
        elif (type(value[0]) == list):
            if (type(value[0][0]) == int):
                return f"extern vector<vector<int>> {key};"
            elif (type(value[0][0]) == float):
                return f"extern vector<vector<float>> {key};"
        else:
            print("Error: Unsupported type in config file (list section)")
            exit(1)

def format_func_line(key, value, section):
    if (type(value) == list):
        if (type(value[0]) == list):
            return f"    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {key}[i][j] = c_{key}[i][j];"
        return f"    for (int i = 0; i < 3; i++) {key}[i] = c_{key}[i];"
    elif (section == 'BANDS'):
        return f"    for (int i = 0; i < nbnd; i++) {key}.push_back(c_{key}[i]);"
    return f"    {key} = c_{key};"

start_phrase = '# Start variable definitions'
end_phrase = '# End variable definitions'

start_func_phrase = '    // Load the C++ configuration file'
end_func_phrase = '    // End of Global Functions'

start_config_phrase = '    // Global Variables in Config struct'
end_config_phrase = '    // End of Global Config Variables'

def write_julia(ALL):
    file_path = 'load/julia_config.jl'
    add_lines(ALL, file_path, start_phrase, end_phrase, format_var_line)
    #add_lines(ALL, file_path, start_func_phrase, end_func_phrase, format_func_line)
    replace_comments(file_path)
    print(f"Successfully updated the file '{file_path}'.")

