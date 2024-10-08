Python API
==========

Building:
---------

Navigate to the src/python directory.

Create and activate a Python virtual environment:
    
    python3 -m venv venv
    
    source venv/bin/activate  

or:

    pipenv shell

Ensure pybind11 and numpy are available in your environment. Install it via pip if necessary:
    
    pip install pybind11 numpy

Build and install the Python bindings:
    
    python pyffwdsetup.py build
    
    python pyffwdsetup.py install

Run the Python test script:
    
    python pyffwdtest.py

#TODO Ramishka is this current? ^^^

Dependencies:
-------------

 * pybind11
 * homebrew
 * gsl

Loading:
--------

once installed the module can be loaded from within the virtual environment with

    >>> import pyffwd

Running Unit Tests
------------------
Navigate to the directory
    cd cancer-integrator/src

Grant execution permission to the script
    chmod +x run_unit_tests.sh

Run the unit tests
    ./run_unit_tests.sh

Run on Jupyter
--------------

Change the directory to the python folder
	os.chdir('/path/to/cancer-integrator/src/python')

Install pybind11 
	!pip install pybind11

Build the pyffwdsetup.py
	!python pyffwdsetup.py build_ext --inplace

Import the pyffwd
	import pyffwd


Example (five-stage-characteristic.py) get results in NumPy array
-----------------------------------------------------------------

import numpy as np
import pyffwd 

def main():
    # Create model with 6 stages
    parameters = pyffwd.Model.create_model(
        6, 
        birth_rates=[0, 0, 0.2, 0.27, 0.27, 0],
        death_rates=[0, 0, 0, 0, 0, 0],
        initial_pops=[1e8, 0, 0, 0, 0, 0]
    )

    # System coefficients
    parameters.m_migr = [
        {1: 2.86e-4},
        {2: 1.06e-5},
        {3: 9.00e-7},
        {4: 1.36e-4},
        {5: 4.56e-7}
    ]

    # Initialize qvalues
    default_qvalues = [1, 1, 1, 1, 1, 1]
    qvalues = [pyffwd.RealVector(default_qvalues) for _ in range(parameters.m_stages)]
    
    for vertex in range(parameters.m_stages):
        qvalues[vertex][vertex] = 0

    time = 0.0
    tmax = 80.0
    dt = 1.0

    # Create a list to store the results
    results = []

    while time < tmax:
        subdivision = 20
        dt2 = dt / subdivision
        for vertex in range(parameters.m_stages):
            for _ in range(subdivision):
                pyffwd.heun_q_step(qvalues[vertex], time, dt2, parameters)

        time += dt

        # Compute probabilities and append results to the list
        row = [int(time)]  # Start with time as integer
        for vertex in range(parameters.m_stages):
            prob = pyffwd.generating_function(qvalues[vertex], parameters.m_initial_pops)
            # Limit to 6 significant digits:
            formatted_prob = 1.0 - prob
            row.append(formatted_prob)

        results.append(row)

    # Convert the results to a NumPy array
    results_array = np.array(results)
    return results_array

if __name__ == "__main__":
    output_array = main()
    print(output_array)



Note
----

GSL in Python: GSL (GNU Scientific Library) in Python, ensure that the random number generator is explicitly set to gsl_rng_mt19937 (Mersenne Twister) for consistency with the C++ implementation.

    const gsl_rng_type *T = gsl_rng_mt19937;


Data Type Conversions: Custom conversion functions for types like std::vector<int> and std::vector<real_t> build to avoid unexpected behaviour and ensure reliable data exchange between Python Pybind11 and C++.


Authors and acknowledgements
----------------------------

* Ramishka Dona Liyanage
