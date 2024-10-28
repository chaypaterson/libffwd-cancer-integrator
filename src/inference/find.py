import re

# Regex to match C++ function definitions
func_def_pattern = re.compile(r'^(?!.*\breturn\b)[\w\s*&]+ ([a-zA-Z_][\w]*)\(')
# Regex to match function calls
func_call_pattern = re.compile(r'\b(?:return\s+)?([a-zA-Z_][\w]*)(?=\s*\()')

def parse_cpp_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    defined_functions = set()  # Set to store func_def_pattern
    call_graph = []  # List to store call relationships
    parent_function = None
    inside_function = False

    # Collect all function definitions
    for line in lines:
        line = line.strip()
        func_def_match = func_def_pattern.match(line)
        if func_def_match:
            function_name = func_def_match.group(1)
            defined_functions.add(function_name)

    # Build the call graph using only known functions
    for line in lines:
        line = line.strip()

        # Look for function definitions
        if not inside_function:
            func_def_match = func_def_pattern.match(line)
            if func_def_match:
                parent_function = func_def_match.group(1)
                inside_function = True
                continue
        
        # look for function calls
        if inside_function:
            # Check for function calls, excluding array accesses and undefined functions
            func_calls = func_call_pattern.findall(line)
            for func in func_calls:
                if func in defined_functions and func != parent_function:
                    # Add to call graph as parent -> child relationship
                    call_graph.append(f"{parent_function} -> {func}")

            # Check for closing brace
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
