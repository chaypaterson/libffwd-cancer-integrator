import re

# Regex to match C++ function definitions
func_def_pattern = re.compile(r'^(?!.*\breturn\b)[\w\s*&]+ ([a-zA-Z_][\w]*)\(')

# Regex to match function calls
func_call_pattern = re.compile(r'\b(?:return\s+)?([a-zA-Z_][\w]*)(?=\s*\()')

# Regex to match function pointer use
func_ptr_assign_pattern = re.compile(r'([a-zA-Z_][\w]*)\s*=\s*([a-zA-Z_][\w]*);')

# Regex to match function pointer definitions
func_pointer_def_pattern = re.compile(r'([a-zA-Z_][\w]*)\s*\(\*([a-zA-Z_][\w]*)\)\s*\(.*\);')

def parse_cpp_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    defined_functions = set()
    defined_func_pointers = set()
    call_graph = []
    parent_function = None
    inside_function = False

    # Collect all function definitions and function pointer definitions
    for line in lines:
        line = line.strip()

        # Check for function definitions
        func_def_match = func_def_pattern.match(line)
        if func_def_match:
            function_name = func_def_match.group(1)
            defined_functions.add(function_name)

        # Check for function pointer definitions
        func_ptr_def_match = func_pointer_def_pattern.match(line)
        if func_ptr_def_match:
            pointer_name = func_ptr_def_match.group(2)
            defined_func_pointers.add(pointer_name)

    # Build the call graph using functions and function pointers
    for line in lines:
        line = line.strip()

        # Look for function definitions
        if not inside_function:
            func_def_match = func_def_pattern.match(line)
            if func_def_match:
                parent_function = func_def_match.group(1)
                inside_function = True
                continue

        if inside_function:
            # Check for function calls
            func_calls = func_call_pattern.findall(line)
            for func in func_calls:
                if func in defined_functions and func != parent_function:
                    call_graph.append(f"{parent_function} -> {func}")

            # Check for function pointer
            func_ptr_assign_match = func_ptr_assign_pattern.search(line)
            if func_ptr_assign_match:
                pointer_name = func_ptr_assign_match.group(1)
                assigned_function = func_ptr_assign_match.group(2)
                if assigned_function in defined_functions:
                    call_graph.append(f"{parent_function} -> *{pointer_name} ({assigned_function})")

            # Check for function pointer usage
            for pointer_name in defined_func_pointers:
                if pointer_name in line:
                    call_graph.append(f"{parent_function} -> *{pointer_name}")

            if re.match(r'^\s*}$', line):
                inside_function = False
                parent_function = None

    return call_graph

# Main program
if __name__ == "__main__":
    filename = 'likelihood-optimisation.cpp'
    call_graph = parse_cpp_file(filename)
    
    # Print the result
    for call in call_graph:
        print(call)
