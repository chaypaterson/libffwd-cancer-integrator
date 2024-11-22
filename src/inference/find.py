import re

# Function definitions
func_def_pattern = re.compile(r'^(?!.*\breturn\b)[\w\s*&]+ ([a-zA-Z_][\w]*)\(')
# Match function calls
func_call_pattern = re.compile(r'\b(?:return\s+)?([a-zA-Z_][\w]*)(?=\s*\()')
# Only whitespace and a closing brace
closing_brace_pattern = re.compile(r'^}$')
# Single-line comments
single_line_comment_pattern = re.compile(r'^\s*//')
# Multi-line comments (start)
multi_line_comment_start_pattern = re.compile(r'/\*')
# Multi-line comments (end)
multi_line_comment_end_pattern = re.compile(r'\*/')
# Pointer declarations and assignments
func_pointer_decl_pattern = re.compile(r'\b(?:[\w\s*&]+)\(\*([a-zA-Z_][\w]*)\)\s*\(')
func_pointer_assign_pattern = re.compile(r'\b([a-zA-Z_][\w]*)\s*=\s*([a-zA-Z_][\w]*)\s*;')

def parse_cpp_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    defined_functions = set()
    function_pointers = {}
    parent_function = None
    inside_function = False
    inside_multi_line_comment = False
    result = []

    for line in lines:
        # Ignore lines that are comments
        if single_line_comment_pattern.match(line):
            continue
        if inside_multi_line_comment:
            if multi_line_comment_end_pattern.search(line):
                inside_multi_line_comment = False
            continue
        if multi_line_comment_start_pattern.search(line):
            inside_multi_line_comment = True
            continue

        # function pointer declarations
        func_pointer_decl_match = func_pointer_decl_pattern.search(line)
        if func_pointer_decl_match:
            pointer_name = func_pointer_decl_match.group(1)
            function_pointers[pointer_name] = None
            continue

        # function definitions
        if not inside_function:
            func_def_match = func_def_pattern.match(line)
            if func_def_match:
                func_name = func_def_match.group(1)
                defined_functions.add(func_name)
                parent_function = func_name
                inside_function = True
                continue

        # Inside function body
        if inside_function:
            # Look for function calls
            func_calls = func_call_pattern.findall(line)
            for func in func_calls:
                # Exclude self-calls and undefined function calls
                if parent_function != func and func in defined_functions:
                    result.append(f"{parent_function} -> {func}")

            # Function pointer assignments
            func_pointer_assign_match = func_pointer_assign_pattern.search(line)
            if func_pointer_assign_match:
                pointer_name = func_pointer_assign_match.group(1)
                assigned_function = func_pointer_assign_match.group(2)
                if pointer_name in function_pointers:
                    function_pointers[pointer_name] = assigned_function
                    result.append(f"{parent_function} -> *{pointer_name} ({assigned_function})")

            # Closing brace
            if closing_brace_pattern.match(line):
                inside_function = False
                parent_function = None

    return result

# Main program
if __name__ == "__main__":
    filename = 'likelihood-optimisation.cpp'
    call_graph = parse_cpp_file(filename)

    # Print the result
    for call in call_graph:
        print(call)
