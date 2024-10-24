import re

#match C++ function definitions
func_def_pattern = re.compile(r'^(?!.*\breturn\b)[\w\s*&]+ ([a-zA-Z_][\w]*)\(')
# match function calls
func_call_pattern = re.compile(r'\b(?:return\s+)?([a-zA-Z_][\w]*)(?=\s*\()')
#match lines with only whitespace and a closing brace
closing_brace_pattern = re.compile(r'^}$')

def parse_cpp_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    parent_function = None
    inside_function = False
    result = []

    for line in lines:
        # Look for function definitions
        if not inside_function:
            func_def_match = func_def_pattern.match(line)
            if func_def_match:
                parent_function = func_def_match.group(1)
                inside_function = True
                continue
        
        #inside function body, look for function calls
        if inside_function:
            func_calls = func_call_pattern.findall(line)
            for func in func_calls:
                # Exclude the current function name to avoid recursion
                if parent_function != func:
                    # Add to result as parent -> child relationship
                    result.append(f"{parent_function} -> {func}")

            # If we encounter a closing brace
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
