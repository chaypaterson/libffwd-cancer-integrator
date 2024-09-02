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

Authors and acknowledgements
----------------------------

* Ramishka Dona Liyanage
