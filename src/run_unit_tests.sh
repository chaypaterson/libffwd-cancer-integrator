#!/bin/bash

# Paths
CORE_DIR="core"
INCLUDE_DIR="/Users/user/cancer-integrator/include"  # Absolute path to the include directory
CPP_TEST_FILES=("tests/five-stage-characteristic.cpp" "tests/likelihood-unit-test.cpp" \
    "tests/ts-loss-characteristic.cpp" "tests/two-hit-characteristic.cpp")
PYTHON_TEST_FILES=("python/Test/five-stage-characteristic.py" "python/Test/likelihood-unit-test.py" \
    "python/Test/ts-loss-characteristic.py" "python/Test/two-hit-characteristic.py")

# Arguments for specific tests
TEST_ARGS=("" "" "0.01 3" "1 10")

# Function to compile and run a C++ test
run_cpp_test() {
    local cpp_test_path=$1
    local test_name=$(basename "$cpp_test_path" .cpp)
    local output_file="cpp_output_${test_name}.txt"
    local args=$2

    echo "Compiling C++ test: $test_name..."
    g++ -std=c++11 -o "$test_name" "$cpp_test_path" \
        -I "$INCLUDE_DIR" -L "$CORE_DIR" -lfast_forward

    if [ $? -ne 0 ]; then
        echo "C++ compilation for $test_name failed!"
        exit 1
    fi

    echo "Running C++ test: $test_name with arguments: $args..."
    ./"$test_name" $args > "$output_file"

    if [ $? -ne 0 ]; then
        echo "C++ test $test_name failed!"
        exit 1
    fi

    echo "C++ test $test_name completed. Output saved to $output_file."
}

# Function to run a Python test
run_python_test() {
    local python_test_path=$1
    local test_name=$(basename "$python_test_path" .py)
    local output_file="python_output_${test_name}.txt"
    local args=$2

    echo "Running Python test: $test_name with arguments: $args..."
    python3 "$python_test_path" $args > "$output_file"

    if [ $? -ne 0 ]; then
        echo "Python test $test_name failed!"
        exit 1
    fi

    echo "Python test $test_name completed. Output saved to $output_file."
}

# Function to compare C++ and Python outputs
compare_outputs() {
    local cpp_output=$1
    local python_output=$2
    local test_name=$(basename "$cpp_output" .txt | sed 's/cpp_output_//')

    echo "Comparing outputs for test: $test_name..."
    diff -u "$cpp_output" "$python_output" > "diff_output_${test_name}.txt"

    if [ $? -eq 0 ]; then
        echo "***Outputs for $test_name are identical***"
        echo "   "
    else
        echo "***Outputs for $test_name differ. See diff_output_${test_name}.txt***"
        echo "   "
    fi
}

# Check if core directory exists
if [ ! -d "$CORE_DIR" ]; then
    echo "Core directory not found!"
    exit 1
fi

# Compile the source files into object files
echo "Compiling fast-forward.cpp..."
cd "$CORE_DIR" || exit
if [ -f "fast-forward.cpp" ]; then
    g++ -std=c++11 -c fast-forward.cpp -I "$INCLUDE_DIR" -o fast-forward.o
else
    echo "Source file fast-forward.cpp not found!"
    exit 1
fi

# Create a static library from the object files
echo "Creating static library libfast_forward.a..."
ar rcs libfast_forward.a fast-forward.o
cd ..

# Run the C++ and Python tests, and compare outputs
for i in "${!CPP_TEST_FILES[@]}"; do
    cpp_test="${CPP_TEST_FILES[$i]}"
    python_test="${PYTHON_TEST_FILES[$i]}"

    # Run the tests with their respective arguments
    run_cpp_test "$cpp_test" "${TEST_ARGS[$i]}"
    run_python_test "$python_test" "${TEST_ARGS[$i]}"

    # Compare outputs
    test_name=$(basename "$cpp_test" .cpp)
    compare_outputs "cpp_output_${test_name}.txt" "python_output_${test_name}.txt"
done
