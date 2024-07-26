from setuptools import setup, Extension
import pybind11

# Path to GSL headers and library
gsl_include_dir = '/opt/homebrew/include'
gsl_lib_dir = '/opt/homebrew/lib'

# Define the extension module
extension_mod = Extension(
    'pyffwd',
    sources=[
        'pyffwd.cpp',  
        '../../src/core/fast-forward.cpp',  
        '../../src/core/graph-model-spec.cpp',
        '../../src/core/gillespie-algorithm.cpp',
    ],
    include_dirs=[
        pybind11.get_include(),
        pybind11.get_include(user=True),
        '../../include',  # Adjust this path
        gsl_include_dir,
        '/opt/homebrew/opt/python@3.12/Frameworks/Python.framework/Versions/3.12/include/python3.12'
    ],
    library_dirs=[
        gsl_lib_dir
    ],
    libraries=['gsl', 'gslcblas'],

    language='c++',
    extra_compile_args=['-std=c++17']
)

# Setup configuration
setup(
    name='pyffwd',
    version='1.0',
    ext_modules=[extension_mod],
    install_requires=['pybind11'],
)

