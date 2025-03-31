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

