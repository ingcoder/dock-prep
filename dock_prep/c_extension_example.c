/*
 * C EXTENSION TEMPLATE FOR PYTHON OPTIMIZATION
 * --------------------------------------------
 * 
 * This file serves as a template/example for creating C extensions to optimize performance-critical
 * parts of Python code. Python is an excellent language for rapid development, but certain
 * computationally intensive operations can benefit significantly from being implemented in C.
 * 
 * When to use C extensions:
 * - For CPU-bound operations where Python's interpreter overhead is a bottleneck
 * - For tight loops or numerical calculations that need maximum performance
 * - When interfacing with existing C libraries
 * 
 * Potential speedup:
 * - Well-designed C extensions can be 10-100x faster than equivalent Python code
 * - Most beneficial for algorithmic or numerical processing
 * 
 * How to use this template:
 * 1. Replace the example functions with your optimized C implementations
 * 2. Update the PyMethodDef array with your function definitions
 * 3. Modify the module name and documentation in PyModuleDef
 * 4. Update the module initialization function name to match your module
 * 5. Compile with appropriate setup.py configuration
 * 
 * For more information on Python C API:
 * https://docs.python.org/3/c-api/
 */

#include <Python.h>
// A simple Python extension module that prints "Hello, World!"

static PyObject* hello_world(PyObject* self, PyObject* noargs) {
    return PyUnicode_FromString("Hello World");
}

// Module methods
static PyMethodDef HelloMethods[] = {
    {"hello_world", (PyCFunction)hello_world, METH_NOARGS, "Prints 'Hello, World!'"},
    {NULL, NULL, 0, NULL}
};

// Module definition
static PyModuleDef hellomodule = {
    PyModuleDef_HEAD_INIT,
    "c_tutorial", // module name
    "A simple module that says hello", // module documentation
    -1, // size of per-interpreter state of the module, or -1 if the module keeps state in global variables
    HelloMethods // methods defined in this module
};

// Module initialization function
PyMODINIT_FUNC PyInit_c_tutorial(void) {
    return PyModule_Create(&hellomodule); // takes the address to the module definition and creates the module
}

