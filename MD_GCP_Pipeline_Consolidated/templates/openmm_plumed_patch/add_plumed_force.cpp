#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "PlumedForce.h"
#include "openmm/System.h"
#include "openmm/Force.h"

// Function to add PlumedForce to System
// This bypasses SWIG entirely by directly calling C++ functions
static PyObject* add_plumed_force_to_system(PyObject* self, PyObject* args) {
    PyObject* system_obj;
    PyObject* plumed_force_obj;
    
    if (!PyArg_ParseTuple(args, "OO", &system_obj, &plumed_force_obj)) {
        return NULL;
    }
    
    // Get the underlying C++ pointers from SWIG objects
    // Both should be SwigPyObject instances
    // We need to access the 'this' attribute which contains the SwigPyObject
    
    // Try to get the 'this' attribute from both objects
    PyObject* system_this = PyObject_GetAttrString(system_obj, "this");
    PyObject* plumed_force_this = PyObject_GetAttrString(plumed_force_obj, "this");
    
    if (!system_this || !plumed_force_this) {
        Py_XDECREF(system_this);
        Py_XDECREF(plumed_force_this);
        PyErr_SetString(PyExc_TypeError, "Objects must have 'this' attribute (SWIG objects)");
        return NULL;
    }
    
    // The 'this' attribute is a SwigPyObject
    // We need to extract the actual C++ pointer from it
    // SwigPyObject has a 'ptr' member that contains the C++ pointer
    
    // Unfortunately, we can't directly access SwigPyObject structure from Python C API
    // We need to use SWIG's runtime functions
    
    // Alternative approach: Use SWIG's Python C API functions
    // But we don't have access to SWIG's runtime in a standalone extension
    
    // Better approach: Use ctypes to call the C++ function directly
    // Or use the fact that SWIG objects can be converted to long (pointer value)
    
    // SWIG objects can be converted to Python int/long which gives the pointer value
    // Use PyNumber_Long to convert, then PyLong_AsVoidPtr to get the pointer
    PyObject* system_ptr_long = PyNumber_Long(system_this);
    PyObject* plumed_force_ptr_long = PyNumber_Long(plumed_force_this);
    
    Py_DECREF(system_this);
    Py_DECREF(plumed_force_this);
    
    if (!system_ptr_long || !plumed_force_ptr_long || PyErr_Occurred()) {
        Py_XDECREF(system_ptr_long);
        Py_XDECREF(plumed_force_ptr_long);
        PyErr_SetString(PyExc_TypeError, "Could not convert SWIG objects to pointer values");
        return NULL;
    }
    
    // Get the pointer values - PyLong_AsVoidPtr converts Python long to void*
    void* system_ptr_val = PyLong_AsVoidPtr(system_ptr_long);
    void* plumed_force_ptr_val = PyLong_AsVoidPtr(plumed_force_ptr_long);
    
    Py_DECREF(system_ptr_long);
    Py_DECREF(plumed_force_ptr_long);
    
    if (PyErr_Occurred()) {
        PyErr_SetString(PyExc_ValueError, "Error extracting pointer values");
        return NULL;
    }
    
    if (!system_ptr_val || !plumed_force_ptr_val) {
        PyErr_SetString(PyExc_ValueError, "Invalid pointer values (NULL)");
        return NULL;
    }
    
    // Cast the pointer values back to C++ pointers
    OpenMM::System* system = reinterpret_cast<OpenMM::System*>(system_ptr_val);
    PlumedPlugin::PlumedForce* plumed_force = reinterpret_cast<PlumedPlugin::PlumedForce*>(plumed_force_ptr_val);
    
    if (!system || !plumed_force) {
        PyErr_SetString(PyExc_ValueError, "Invalid pointers");
        return NULL;
    }
    
    // Now call the C++ function directly
    try {
        int index = system->addForce(static_cast<OpenMM::Force*>(plumed_force));
        return PyLong_FromLong(index);
    } catch (...) {
        PyErr_SetString(PyExc_RuntimeError, "Exception occurred in C++ addForce");
        return NULL;
    }
}

// Method definitions
static PyMethodDef module_methods[] = {
    {
        "add_plumed_force_to_system",
        add_plumed_force_to_system,
        METH_VARARGS,
        "Add a PlumedForce to an OpenMM System (bypasses SWIG)"
    },
    {NULL, NULL, 0, NULL}  // Sentinel
};

// Module definition
static struct PyModuleDef module_def = {
    PyModuleDef_HEAD_INIT,
    "add_plumed_force",  // Module name
    "C extension to add PlumedForce to System",  // Module docstring
    -1,
    module_methods
};

// Module initialization
PyMODINIT_FUNC PyInit_add_plumed_force(void) {
    return PyModule_Create(&module_def);
}

