# exclude_functions.py

# List of C++ keywords and common function names to exclude
excluded_keywords = {
    "alignas", "alignof", "and", "and_eq", "asm", "atomic_cancel", "atomic_commit", 
    "atomic_noexcept", "auto", "bitand", "bitor", "bool", "break", "case", "catch", 
    "char", "char8_t", "char16_t", "char32_t", "class", "compl", "concept", "const", 
    "consteval", "constexpr", "constinit", "const_cast", "continue", "co_await", 
    "co_return", "co_yield", "decltype", "default", "delete", "do", "double", 
    "dynamic_cast", "else", "inverse", "enum", "explicit", "export", "extern", "false", "float", 
    "for", "friend", "goto", "if", "inline", "int", "long", "mutable", "namespace", 
    "new", "noexcept", "not", "not_eq", "nullptr", "operator", "or", "or_eq", "private", 
    "protected", "public", "reflexpr", "register", "reinterpret_cast", "requires", 
    "return", "short", "signed", "sizeof", "static", "static_assert", "static_cast", 
    "struct", "switch", "synchronized", "template", "this", "thread_local", "throw", 
    "true", "try", "typedef", "typeid", "typename", "union", "unsigned", "using", 
    "virtual", "void", "volatile", "wchar_t", "while", "xor", "xor_eq", "type", "real_t",
    "objective"
}

# Dictionary of excluded functions based on library
excluded_functions_by_lib = {
    "cstdio": {
        # Operations on files
        "remove", "rename", "tmpfile", "tmpnam",
        
        # File access
        "fclose", "fflush", "fopen", "freopen", "setbuf", "setvbuf",
        
        # Formatted input/output
        "fprintf", "fscanf", "printf", "scanf", "snprintf", "sprintf",
        "sscanf", "vfprintf", "vfscanf", "vprintf", "vscanf", "vsnprintf", 
        "vsprintf", "vsscanf",
        
        # Character input/output
        "fgetc", "fgets", "fputc", "fputs", "getc", "getchar", "gets", 
        "putc", "putchar", "puts", "ungetc",
        
        # Direct input/output
        "fread", "fwrite",
        
        # File positioning
        "fgetpos", "fseek", "fsetpos", "ftell", "rewind",
        
        # Error-handling
        "clearerr", "feof", "ferror", "perror",
        
        # Macros
        "BUFSIZ", "EOF", "FILENAME_MAX", "FOPEN_MAX", "L_tmpnam", "NULL", "TMP_MAX",
        "_IOFBF", "_IOLBF", "_IONBF", "SEEK_CUR", "SEEK_END", "SEEK_SET",
        
        # Types
        "FILE", "fpos_t", "size_t"
    },

    "cstdlib": {
        # Types
        "div_t", "ldiv_t", "lldiv_t", "size_t",
        
        # Macro constants
        "EXIT_SUCCESS", "EXIT_FAILURE", "MB_CUR_MAX", "NULL", "RAND_MAX",
        
        # Functions - Process control
        "abort", "exit", "quick_exit", "_Exit", "atexit", "at_quick_exit", 
        "system", "getenv",
        
        # Functions - Memory management
        "malloc", "aligned_alloc", "calloc", "realloc", "free",
        
        # Functions - Numeric string conversion
        "atof", "atoi", "atol", "atoll", "strtol", "strtoll", 
        "strtoul", "strtoull", "strtof", "strtod", "strtold",
        
        # Functions - Wide string manipulation
        "mblen", "mbtowc", "wctomb", "mbstowcs", "wcstombs",
        
        # Functions - Miscellaneous algorithms and math
        "rand", "srand", "qsort", "bsearch", 
        "abs", "labs", "llabs", "div", "ldiv", "lldiv"
    },

    "iostream": {
        "std::cin", "std::cout", "std::cerr", "std::clog", "std::getline", 
        "std::flush", "std::setw", "std::setprecision", "std::fixed", 
        "std::scientific", "std::endl", "std::hex", "std::dec", "std::oct"
    },

    "fstream": {
        # Types
        "std::fstream", "std::wfstream",
        
        # Member types
        "char_type", "traits_type", "int_type", "pos_type", "off_type", 
        "native_handle_type(C++26)",
        
        # Member functions
        "basic_fstream", "~basic_fstream", "operator=", "swap",
        "rdbuf", "native_handle(C++26)",
        
        # File operations
        "is_open", "open", "close",
        
        # Inherited functions
        "operator>>", "get", "peek", "unget", "putback", "getline", 
        "ignore", "read", "readsome", "gcount", 
        "tellg", "seekg", "sync",
        "operator<<", "put", "write", "tellp", "seekp", "flush",
        "good", "eof", "fail", "bad", "operator!", "operator bool", 
        "rdstate", "setstate", "clear", "copyfmt", "fill", 
        "exceptions", "imbue", "rdbuf", "tie", "narrow", "widen",
        "flags", "setf", "unsetf", "precision", "width", 
        "imbue", "getloc", "xalloc", "iword", "pword", 
        "register_callback", "sync_with_stdio(static)",
        
        # Non-member functions
        "std::swap(std::basic_fstream)",
    },

    "vector": {
        # Member functions
        "push_back", "pop_back", "size", "empty", "clear", "resize", 
        "capacity", "at", "operator[]", "front", "back", "reserve", 
        "shrink_to_fit", "insert", "erase", "emplace", "emplace_back", 
        "swap", "assign", "get_allocator", "data",
        
        # Iterators
        "begin", "end", "rbegin", "rend", "cbegin(C++11)", "cend(C++11)", 
        "crbegin(C++11)", "crend(C++11)",
        
        # Element Access
        "at", "operator[]", "front", "back", "data",
    
        # Capacity functions
        "empty", "size", "max_size", "reserve", "capacity", "shrink_to_fit",
    
        # Modifiers
        "clear", "insert", "insert_range(C++23)", "emplace", "erase", 
        "push_back", "emplace_back", "append_range(C++23)", "pop_back", 
        "resize", "swap"
    },

    "cmath": {
        "abs(float)", "fabs", "fabsf", "fabsl", "fmod", "fmodf", "fmodl",
        "remainder", "remainderf", "remainderl", "remquo", "remquof", "remquol",
        "fma", "fmaf", "fmal", "fmax", "fmaxf", "fmaxl", "fmin", "fminf", "fminl",
        "fdim", "fdimf", "fdiml", "nan", "nanf", "nanl", "lerp(C++20)",
        "exp", "expf", "expl", "exp2", "exp2f", "exp2l", "expm1", "expm1f", "expm1l",
        "log", "logf", "logl", "log10", "log10f", "log10l", "log2", "log2f", "log2l",
        "log1p", "log1pf", "log1pl", "pow", "powf", "powl", "sqrt", "sqrtf", "sqrtl",
        "cbrt", "cbrtf", "cbrtl", "hypot", "hypotf", "hypotl", "sin", "sinf", "sinl",
        "cos", "cosf", "cosl", "tan", "tanf", "tanl", "asin", "asinf", "asinl", 
        "acos", "acosf", "acosl", "atan", "atanf", "atanl", "atan2", "atan2f", 
        "atan2l", "sinh", "sinhf", "sinhl", "cosh", "coshf", "coshl", "tanh", "tanhf", 
        "tanhl", "asinh", "asinhf", "asinhl", "acosh", "acoshf", "acoshl", "atanh", 
        "atanhf", "atanhl", "erf", "erff", "erfl", "erfc", "erfcf", "erfcl", "lgamma", 
        "lgammaf", "lgammal", "tgamma", "tgammaf", "tgammal", "ceil", "ceilf", "ceill",
        "floor", "floorf", "floorl", "trunc", "truncf", "truncl", "round", "roundf", 
        "roundl", "rint", "rintf", "rintl", "nearbyint", "nearbyintf", "nearbyintl", 
        "remquo", "remquof", "remquol", "modf", "modff", "modfl", "scalbn", "scalbnf", 
        "scalbnl", "scalbln", "scalblnf", "scalblnl", "logb", "logbf", "logbl", 
        "frexp", "frexpf", "frexpl", "ldexp", "ldexpf", "ldexpl", "ilogb", "ilogbf", 
        "ilogbl", "copysign", "copysignf", "copysignl", "nan", "nanf", "nanl",
        "signbit", "signbitf", "signbitl", "isnan", "isnanf", "isnanl", "isfinite", 
        "isfinitef", "isfinitel", "isinf", "isinff", "isinfl", "isnormal", "isnormalf", 
        "isnormall", "isgreater", "isgreaterf", "isgreaterl", "isless", "islessf", 
        "islessl", "islessequal", "islessequalf", "islessequall", "islessgreater", 
        "islessgreaterf", "islessgreaterl", "isunordered", "isunorderedf", "isunorderedl"
    },
    
    "functional": {
        "std::function", "std::bind", "std::placeholders", "std::less", 
        "std::greater", "std::equal_to", "std::negate"
    },
    
    "string": {
        "append", "size", "length", "empty", "substr", "find", 
        "replace", "erase", "compare", "c_str", "operator[]", "at", 
        "front", "back", "data", "operator basic_string_view", 
        "begin", "end", "rbegin", "rend", "reserve", "capacity", 
        "shrink_to_fit", "clear", "insert", "push_back", "pop_back", 
        "operator+", "operator==", "operator!=", "operator<", "operator>", 
        "operator<=", "operator>=", "std::swap", "getline", 
        "stoi", "stol", "stoll", "stoul", "stoull", "stof", "stod", 
        "stold", "to_string", "to_wstring", "insert_range", "append_range", 
        "replace_with_range", "resize", "resize_and_overwrite", "find_first_of", 
        "find_first_not_of", "find_last_of", "find_last_not_of", 
        "starts_with", "ends_with", "contains", "erase_if", "memcpy"
        "memcpy_s", "std::memcpy"
    },

    "thread": {
        "std::thread", "std::mutex", "std::lock_guard", "std::unique_lock", 
        "std::async", "std::this_thread::sleep_for", "std::this_thread::sleep_until"
        "join", "thread"
    }

}

# Function to return all excluded functions as a set
def get_excluded_functions():
    excluded = set(excluded_keywords)  # Add the general keywords first
    # Add the functions from the libraries
    for lib_functions in excluded_functions_by_lib.values():
        excluded.update(lib_functions)
    return excluded
