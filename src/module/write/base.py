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


def remove_multiple_empty_lines(file_path):
    with open(file_path, "r") as file:
        lines = file.readlines()
    updated_lines = []
    for i in range(len(lines)):
        if i < len(lines) - 1:
            if lines[i] == "\n" and lines[i + 1] == "\n":
                continue
        updated_lines.append(lines[i])
    with open(file_path, "w") as file:
        file.writelines(updated_lines)
