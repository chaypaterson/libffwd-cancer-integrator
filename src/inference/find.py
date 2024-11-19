import re
from exclude_functions import get_excluded_functions  # Import the exclusion list

# Match C++ function definitions
func_def_pattern = re.compile(r'^(?!.*\breturn\b)[\w\s*&]+ ([a-zA-Z_][\w]*)\(')
# Match function calls
func_call_pattern = re.compile(r'\b(?:return\s+)?([a-zA-Z_][\w]*)\s*\(')
# Match lines with only whitespace and a closing brace
closing_brace_pattern = re.compile(r'^}$')
# Match single-line comments
single_line_comment_pattern = re.compile(r'^\s*//')
# Match multi-line comments (start)
multi_line_comment_start_pattern = re.compile(r'/\*')
# Match multi-line comments (end)
multi_line_comment_end_pattern = re.compile(r'\*/')

def parse_cpp_file(filename, excluded_functions):
    with open(filename, 'r') as file:
        lines = file.readlines()

    parent_function = None
    inside_function = False
    inside_multi_line_comment = False
    result = []

    for line in lines:
        # Ignore lines that are comments or are inside a comment block
        if single_line_comment_pattern.match(line):
            continue
        if inside_multi_line_comment:
            if multi_line_comment_end_pattern.search(line):
                inside_multi_line_comment = False
            continue
        if multi_line_comment_start_pattern.search(line):
            inside_multi_line_comment = True
            continue
        
        # Look for function definitions (not inside a comment)
        if not inside_function:
            func_def_match = func_def_pattern.match(line)
            if func_def_match:
                parent_function = func_def_match.group(1)
                inside_function = True
                continue
        
        # Inside function body, look for function calls
        if inside_function:
            func_calls = func_call_pattern.findall(line)
            for func in func_calls:
                # Exclude functions from known libraries and the parent function
                if func not in excluded_functions and parent_function != func:
                    result.append(f"{parent_function} -> {func}")

            # If we encounter a closing brace, end the current function
            if closing_brace_pattern.match(line):
                inside_function = False
                parent_function = None

    return result

# Main program
if __name__ == "__main__":
    # Get the set of excluded functions from standard libraries
    excluded_functions = get_excluded_functions()

    filename = 'likelihood-optimisation.cpp'
    call_graph = parse_cpp_file(filename, excluded_functions)
    
    # Print the result
    for call in call_graph:
        print(call)
