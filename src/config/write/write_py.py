from .write_c import add_lines

start_phrase = "### Variables ###"
end_phrase = "### End Variables ###"

start_set_phrase = "            # Set the variable"
end_set_phrase = "            # Finished setting variables"

start_init_phrase = "# Begin python->c++ interface"
end_init_phrase = "# End python->c++ interface"

start_module_phrase = "    // Begin the Config class"
end_module_phrase = "    // End the Config class"

def format_var_line(key, value, section):
    if type(value) == str:
        value = "'" + value + "'"
    if section == 'BANDS':
        return f"{key} = []\n{key}.append({value})"
    return f"{key} = {value}"

def format_module_line(key, value, section):
    return f"    .def_readwrite(\"{key}\", &Config::{key})"

def format_init_line(key, value, section):
    return f"placeholder.{key} = config.{key}"

def format_set_line(key, value, section):
    if type(value) == str:
        if section == 'BANDS':
            return f"            if \"{key}\" in key:\n                global {key}\n                {key}.append(value)"
        return f"            if \"{key}\" in key:\n                global {key}\n                {key} = value"
    if type(value) == int:
        if section == 'BANDS':
            return f"            if \"{key}\" in key:\n                global {key}\n                {key}.append(int(value))"
        if key == 'dimension':
            return f"            if \"{key}\" in key:\n                global {key}\n                {key} = int(value)\n                got_dimension = True"
        return f"            if \"{key}\" in key:\n                global {key}\n                {key} = int(value)"
    if type(value) == float:
        if section == 'BANDS':
            return f"            if \"{key}\" in key:\n                global {key}\n                {key}.append(float(value))"
        return f"            if \"{key}\" in key:\n                global {key}\n                {key} = float(value)"
    if type(value) == bool:
        if section == 'BANDS':
            return f"            if \"{key}\" in key:\n                global {key}\n                {key}.append(value == 'true')"
        return f"            if \"{key}\" in key:\n                global {key}\n                {key} = value == 'true'"
    if type(value) == list:
        if type(value[0]) == str:
            return f"            if \"{key}\" in key:\n                global {key}\n                {key} = [value.split()[i] for i in range(3)]"
        if type(value[0]) == int:
            return f"            if \"{key}\" in key:\n                global {key}\n                {key} = [int(value.split()[i]) for i in range(3)]"
        if type(value[0]) == float:
            return f"            if \"{key}\" in key:\n                global {key}\n                {key} = [float(value.split()[i]) for i in range(3)]"
        if type(value[0]) == bool:
            return f"            if \"{key}\" in key:\n                global {key}\n                {key} = [value.split()[i] == 'true' for i in range(3)]"
        if type(value[0]) == list:
            return f"            if section == \"{section}\" and index < 3:\n                global {key}\n                {key}.append([float(line.split()[i]) for i in range(3)])\n                index += 1"
    print(f"Error: {key} has an unsupported type.")
    exit()
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

def write_py(ALL):
    file_path = 'load/config.py'
    add_lines(ALL, file_path, start_phrase, end_phrase, format_var_line)
    add_lines(ALL, file_path, start_set_phrase, end_set_phrase, format_set_line)
    replace_comments(file_path)
    print(f"Successfully updated the file '{file_path}'.")
    #file_path = '../fmodule.cpp'
    #add_lines(ALL, file_path, start_module_phrase, end_module_phrase, format_module_line)
    #print(f"Successfully updated the file '{file_path}'.")
    file_path = '../__init__.py'
    add_lines(ALL, file_path, start_init_phrase, end_init_phrase, format_init_line)
    replace_comments(file_path)
