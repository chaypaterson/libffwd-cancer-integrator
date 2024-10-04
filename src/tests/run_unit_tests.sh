#!/bin/bash

# Paths
# TODO currently assumes we are running in src/, maybe change to src/tests?
CORE_DIR="../core"
INCLUDE_DIR="../../include"
GSL_DIR="$(brew --prefix gsl)"
CPP_TEST_FILES=("ts-loss-tau.cpp" "five-stage-characteristic.cpp" "likelihood-unit-test.cpp" \
    "ts-loss-characteristic.cpp" "two-hit-characteristic.cpp" \
    "gillespie-sampler.cpp" "two-hit-gillespie.cpp" \
    "../errors/ts-loss-characteristic-errors.cpp" "ts-loss-gillespie-poly.cpp" \
    "../errors/gillespie-errors.cpp" "ts-loss-gillespie.cpp") 
PYTHON_TEST_FILES=("../python/Test/ts-loss-tau.py" "../python/Test/five-stage-characteristic.py" \
    "../python/Test/likelihood-unit-test.py" "../python/Test/ts-loss-characteristic.py" \
    "../python/Test/two-hit-characteristic.py" \
    "../python/Test/gillespie-sampler.py" "../python/Test/two-hit-gillespie.py" \
    "../python/Test/ts-loss-characteristic-errors.py" "../python/Test/ts-loss-gillespie-poly.py" \
    "../python/Test/gillespie-errors.py" "../python/Test/ts-loss-gillespie.py") 

BINDIR=.bin
LIBDIR=.libs
LIBFFWD=$LIBDIR/libffwd.so
LIBGILL=$LIBDIR/libgillespie.so
OUTDIR=results
mkdir -p $BINDIR
mkdir -p $LIBDIR
mkdir -p $OUTDIR

# Write a note to the summary file:
date >> test_failures.txt
echo "Test failure summary:" >> test_failures.txt 

# Compile the libraries:
echo "Compiling fast-forward.cpp..."
pwd # which directory are we running this in?
ls $CORE_DIR
if [ -f "$CORE_DIR/fast-forward.cpp" ]; then
    g++ -O2 -std=c++11 -c $CORE_DIR/fast-forward.cpp -I "$INCLUDE_DIR" -o $LIBFFWD
else
    echo "Source file fast-forward.cpp not found!"
    exit 1
fi

if [ -f "$CORE_DIR/gillespie-algorithm.cpp" ]; then
    echo "Compiling gillespie-algorithm.cpp..."
    g++ -O2 -std=c++11 -c $CORE_DIR/gillespie-algorithm.cpp -I "$INCLUDE_DIR" -I "$GSL_DIR/include" \
        -o $LIBGILL
else
    echo "Source file gillespie-algorithm.cpp not found!"
    exit 1
fi


# Arguments for specific tests
TEST_ARGS=("1 100 0.1" "" "" "0.01 3" "1 10" "1 100" "1 10" "0.01" "123 10" "1 10" "1 10")

# Function to compile and run a C++ test
run_cpp_test() {
    local cpp_test_path=$1
    local test_name=$(basename "$cpp_test_path" .cpp)
    local output_file="cpp_output_${test_name}.txt"
    local args=$2

    echo "Compiling C++ test: $test_name..."
    g++ -O2 -std=c++11 -o $BINDIR/"$test_name" \
        $LIBFFWD $LIBGILL "$cpp_test_path" \
        -I "$INCLUDE_DIR" -I "$GSL_DIR/include" \
        -L "$CORE_DIR" -L "$GSL_DIR/lib" \
        -lgsl -lgslcblas -pthread

    if [ $? -ne 0 ]; then
        echo "C++ compilation for $test_name failed!"
        exit 1
    fi

    echo "Running C++ test: $test_name with arguments: $args..."
    ./$BINDIR/"$test_name" $args > $OUTDIR/"$output_file"

    if [ $? -ne 0 ]; then
        echo "C++ test $test_name failed!"
        exit 1
    fi

    echo "C++ test $test_name completed. Output saved to $OUTDIR/$output_file."
}

# Function to run a Python test
run_python_test() {
    local python_test_path=$1
    local test_name=$(basename "$python_test_path" .py)
    local output_file="python_output_${test_name}.txt"
    local args=$2

    echo "Running Python test: $test_name with arguments: $args..."
    python3 "$python_test_path" $args > $OUTDIR/"$output_file"

    if [ $? -ne 0 ]; then
        echo "Python test $test_name failed!"
        exit 1
    fi

    echo "Python test $test_name completed. Output saved to $OUTDIR/$output_file."
}

# Function to compare C++ and Python outputs
compare_outputs() {
    local cpp_output=$1
    local python_output=$2
    local test_name=$(basename "$cpp_output" .txt | sed 's/cpp_output_//')

    echo "Comparing outputs for test: $test_name..."
    diff -u "$cpp_output" "$python_output" > "diff_output_${test_name}.txt"

    if [ $? -eq 0 ]; then
        echo -e "\033[32m\U2713\033[0m Outputs for $test_name are identical."
        echo "----------------------------------------------"
        rm -v "diff_output_${test_name}.txt"
    else
        echo -e "\033[37;41;6m\U1F4A5 Outputs for $test_name differ!\033[0m See diff_output_${test_name}.txt"
        echo "-------------------------------------------------------------------------------------"
        echo "diff_output_${test_name}.txt" >> test_failures.txt
    fi
}

# Check if core directory exists
if [ ! -d "$CORE_DIR" ]; then
    echo "Core directory not found!"
    exit 1
fi

# Run the C++ and Python tests, and compare outputs
for i in "${!CPP_TEST_FILES[@]}"; do
    cpp_test="${CPP_TEST_FILES[$i]}"
    python_test="${PYTHON_TEST_FILES[$i]}"

    # Run the tests with arguments
    time run_cpp_test "$cpp_test" "${TEST_ARGS[$i]}"
    time run_python_test "$python_test" "${TEST_ARGS[$i]}"

    # Compare outputs
    test_name=$(basename "$cpp_test" .cpp)
    compare_outputs $OUTDIR/"cpp_output_${test_name}.txt" $OUTDIR/"python_output_${test_name}.txt"
done

# All done: alert the user by ringing the bell
printf "\a"
