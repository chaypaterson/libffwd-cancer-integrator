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

Ensure pybind11 is available in your environment. Install it via pip if necessary:
    
    pip install pybind11

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

Note
----

GSL in Python: GSL (GNU Scientific Library) in Python, ensure that the random number generator is explicitly set to gsl_rng_mt19937 (Mersenne Twister) for consistency with the C++ implementation.

    const gsl_rng_type *T = gsl_rng_mt19937;


Data Type Conversions: Custom conversion functions for types like std::vector<int> and std::vector<real_t> build to avoid unexpected behaviour and ensure reliable data exchange between Python Pybind11 and C++.


Authors and acknowledgements
----------------------------

* Ramishka Dona Liyanage
