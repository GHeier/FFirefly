def remove_lines_between_phrases(lines, start, end):
    updated_lines = []
    add_lines = True
    for line in lines:
        # Skip the lines between the matching phrases
        if start in line:
            add_lines = False
            updated_lines.append(line)
        if add_lines:
            updated_lines.append(line)
        if end in line:
            add_lines = True
            updated_lines.append(line)
    return updated_lines

def add_lines_between_phrases(lines, new_lines, start, end):
    updated_lines = []
    add_lines = True
    for line in lines:
        # Skip the lines between the matching phrases
        if start in line:
            add_lines = False
            updated_lines.append(line)
            updated_lines.extend(new_lines)
        if add_lines:
            updated_lines.append(line)
        if end in line:
            add_lines = True
            updated_lines.append(line)
    return updated_lines

def format_var_line(key, value, section):
    if isinstance(value, list):
        # Check if it's a 2D array
        if all(isinstance(sub, list) for sub in value):
            # Get dimensions and element type for 2D array
            rows = len(value)
            cols = len(value[0]) if rows > 0 else 0
            array_type = (
                "float"
                if any(isinstance(x, float) for row in value for x in row)
                else "int"
            )
            array_elements = ', '.join(
                '{' + ', '.join(map(str, row)) + '}' for row in value
            )
            return f"{array_type} c_{key}[{rows}][{cols}] = {{{array_elements}}};"
        else:
            # Handle 1D array
            size = len(value)
            array_type = "float" if any(isinstance(x, float) for x in value) else "int"
            array_elements = ', '.join(map(str, value))
            return f"{array_type} c_{key}[{size}] = {{{array_elements}}};"
    else:
        if isinstance(value, str):
            # Handle strings
            if section == 'BANDS':
                return (
                        f'char** c_{key} = (char*[]){{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};\n'
                        f'char** get_{key}() {{return (char**)c_{key};}}'
                        )
            return (
                f'char* c_{key} = "{value}";\n'
                f'char* get_{key}() {{return c_{key};}}'
            )
        elif isinstance(value, bool):
            # Handle booleans
            if section == 'BANDS':
                return f"bool c_{key}[50];"
            return f"bool c_{key} = {'true' if value else 'false'};"
        else:
            # Handle other scalar values
            value_type = "float" if isinstance(value, float) else "int"
            if section == 'BANDS':
                return f"{value_type} c_{key}[50];"
            return f"{value_type} c_{key} = {value};"

def format_func_line(key, value, section):
    el = 'else '
    if key == 'category':
        el = ''
    index = ''
    if section == 'BANDS':
        index = '[n]'
    if section == 'CELL' or section == 'BRILLOUIN_ZONE':
        return ''

    if (type(value) == str):
        if (key == "band"):
            return f"            {el}if (strstr(key, \"{key}\") != NULL) {{\n                n = atoi(key + 4)-1;\n                set_string(&c_{key}{index}, value);\n            }}"
        return f"            {el}if (strstr(key, \"{key}\") != NULL) {{\n                set_string(&c_{key}, value);\n            }}"
    elif (type(value) == int):
        if key == "dimension":
            return f"            {el}if (strstr(key, \"{key}\") != NULL) {{\n                c_{key}{index} = atoi(value);\n                 got_dimension = true;\n            }}"
        if key == "nbnd":
            return f"            {el}if (strstr(key, \"{key}\") != NULL) {{\n                c_{key}{index} = atoi(value);\n                 got_nbnd = true;\n            }}"
        return f"            {el}if (strstr(key, \"{key}\") != NULL) {{\n                c_{key}{index} = atoi(value);\n            }}"
    elif (type(value) == float):
        return f"            {el}if (strstr(key, \"{key}\") != NULL) {{\n                c_{key}{index} = atof(value);\n            }}"
    elif (type(value) == bool):
        return f"            {el}if (strstr(key, \"{key}\") != NULL) {{\n                strip_single_quotes(value);\n                if (strcmp(value, \"true\") == 0) {{\n                    c_{key}{index} = true;\n                }} else {{\n                    c_{key}{index} = false;\n                }}\n            }}"
    elif (type(value) == list):
        if (type(value[0]) == int):
            return f"            {el}if (strstr(key, \"{key}\") != NULL) {{\n                sscanf(line, \" {key} = %d %d %d\", &c_{key}[0], &c_{key}[1], &c_{key}[2]);\n            }}"
        else:
            print("Key, value, section: ", key, value, section)
            print("Error: Unsupported type in config file (list section)")
            exit(1)
    else:
        print("Error: Unsupported type in config file")
        exit(1)

def format_header_line(key, value, section):
    index = ''
    if section == 'BANDS':
        index = '[50]'
    if (type(value) == str):
        if (key == "band"):
            return f"extern char** c_{key}; char** get_{key}();"
        return f"extern char* c_{key}; char* get_{key}();"
    elif (type(value) == int):
        return f"extern int c_{key}{index};"
    elif (type(value) == float):
        return f"extern float c_{key}{index};"
    elif (type(value) == bool):
        return f"extern bool c_{key}{index};"
    elif (type(value) == list):
        if (type(value[0]) == int):
            return f"extern int c_{key}{index}[3];"
        elif (type(value[0]) == float):
            return f"extern float c_{key}{index}[3];"
        elif (type(value[0]) == list):
            if (type(value[0][0]) == int):
                return f"extern int c_{key}{index}[3][3];"
            elif (type(value[0][0]) == float):
                return f"extern float c_{key}{index}[3][3];"
        else:
            print("Error: Unsupported type in config file (list section)")
            exit(1)
    else:
        print("Error: Unsupported type in config file")
        exit(1)

def format_unload_line(key, value, section):
    if type(value) == str:
        return f'    free(c_{key});'
    if key == 'band':
        return f'    for (int i = 0; i < 50; i++) {{\n        free(c_{key}[i]);\n    }}'
    return ''


def get_new_lines(ALL, func):
    new_lines = []
    for section in ALL:
        new_lines.append('\n//' + '[' + section + ']\n')
        for key, value in ALL[section].items():
            new_lines.append(func(key, value, section) + '\n')
        new_lines.append('')
    return new_lines

def add_lines(ALL, file_path, start_phrase, end_phrase, func):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    lines = remove_lines_between_phrases(lines, start_phrase, end_phrase)
    new_lines = get_new_lines(ALL, func)
    updated_lines = add_lines_between_phrases(lines, new_lines, start_phrase, end_phrase)

    # Write the modified content back to the file
    with open(file_path, 'w') as file:
        file.writelines(updated_lines)

def remove_multiple_empty_lines(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    updated_lines = []
    for i in range(len(lines)):
        if i < len(lines) - 1:
            if lines[i] == '\n' and lines[i + 1] == '\n':
                continue
        updated_lines.append(lines[i])
    with open(file_path, 'w') as file:
        file.writelines(updated_lines)

start_phrase = '// Global Variables are listed below, with their default values'
end_phrase = '// End of Global Variables'

start_func_phrase = '            // Read in variable values from the config file'
end_func_phrase = '            // End of variable reading'

start_unload_phrase = '// Unload the config file'
end_unload_phrase = '// End of unloading the config file'

def write_c_header(ALL):
    file_path = 'load/c_config.h'
    add_lines(ALL, file_path, start_phrase, end_phrase, format_header_line)
    remove_multiple_empty_lines(file_path)
    print(f"Successfully updated the file '{file_path}'.")

def write_c(ALL):
    file_path = 'load/c_config.c'
    add_lines(ALL, file_path, start_phrase, end_phrase, format_var_line)
    add_lines(ALL, file_path, start_func_phrase, end_func_phrase, format_func_line)
    add_lines(ALL, file_path, start_unload_phrase, end_unload_phrase, format_unload_line)
    remove_multiple_empty_lines(file_path)
    print(f"Successfully updated the file '{file_path}'.")
