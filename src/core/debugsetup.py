from setuptools import setup, Extension
import pybind11

# Define the extension module
extension_mod = Extension(
    'pybinding',
    sources=['debugbindding.cpp', 'fast-forward.cpp', 'graph-model-spec.cpp'],
    include_dirs=[
        pybind11.get_include(),
        pybind11.get_include(user=True),
        '/Users/user/cancer-integrator/include',  # Adjust this path
        '/opt/homebrew/opt/python@3.12/Frameworks/Python.framework/Versions/3.12/include/python3.12'
    ],
    language='c++',
    extra_compile_args=['-std=c++17', '-stdlib=libc++', '-mmacosx-version-min=10.9']
)

# Setup configuration
setup(
    name='pybinding',
    version='1.0',
    ext_modules=[extension_mod],
    install_requires=['pybind11'],
)

