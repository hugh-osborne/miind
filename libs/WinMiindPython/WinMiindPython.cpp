#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <WinMiindLib/SimulationParserCPU.h>

SimulationParserCPU<MPILib::CustomConnectionParameters>* model;

PyObject* miind_init(PyObject* self, PyObject* args)
{
    if (model) {
        delete model;
        model = NULL;
    }

    int nodes;
    char* filename;
    if (PyArg_ParseTuple(args, "is", &nodes, &filename))
        model = new SimulationParserCPU<MPILib::CustomConnectionParameters>(nodes, std::string(filename));
    else if (PyArg_ParseTuple(args, "s", &filename))
        model = new SimulationParserCPU<MPILib::CustomConnectionParameters>(std::string(filename));
    else
        return NULL;
        
    model->init();

    Py_RETURN_NONE;
}

PyObject* miind_getTimeStep(PyObject* self, PyObject* args)
{
    return Py_BuildValue("d", model->getTimeStep());
}

PyObject* miind_getSimulationLength(PyObject* self, PyObject* args)
{
    return Py_BuildValue("d", model->getSimulationLength());
}

PyObject* miind_startSimulation(PyObject* self, PyObject* args)
{
    model->startSimulation();

    Py_RETURN_NONE;
}

PyObject* miind_evolveSingleStep(PyObject* self, PyObject* args)
{
    PyObject* float_list;
    int pr_length;

    if (!PyArg_ParseTuple(args, "O", &float_list))
        return NULL;
    pr_length = PyObject_Length(float_list);
    if (pr_length < 0)
        return NULL;

    std::vector<double> activities(pr_length);

    for (int index = 0; index < pr_length; index++) {
        PyObject* item;
        item = PyList_GetItem(float_list, index);
        if (!PyFloat_Check(item))
            activities[index] = 0.0;
        activities[index] = PyFloat_AsDouble(item);
    }

    std::vector<double> out_activities = model->evolveSingleStep(activities);

    PyObject* tuple = PyTuple_New(out_activities.size());

    for (int index = 0; index < out_activities.size(); index++) {
        PyTuple_SetItem(tuple, index, Py_BuildValue("d", out_activities[index]));
    }

    return tuple;
}

PyObject* miind_endSimulation(PyObject* self, PyObject* args)
{
    model->endSimulation();

    Py_RETURN_NONE;
}

/*
 * List of functions to add to WinMiindPython in exec_WinMiindPython().
 */
static PyMethodDef WinMiindPython_functions[] = {
    {"init",  miind_init, METH_VARARGS, "Init Miind Model."},
    {"getTimeStep",  miind_getTimeStep, METH_VARARGS, "Get time step."},
    {"getSimulationLength",  miind_getSimulationLength, METH_VARARGS, "Get sim time."},
    {"startSimulation",  miind_startSimulation, METH_VARARGS, "Start simulation."},
    {"evolveSingleStep",  miind_evolveSingleStep, METH_VARARGS, "Evolve one time step."},
    {"endSimulation",  miind_endSimulation, METH_VARARGS, "Clean up."},
    { NULL, NULL, 0, NULL } /* marks end of array */
};

/*
 * Initialize WinMiindPython. May be called multiple times, so avoid
 * using static state. "oops." - Hugh
 */
int exec_WinMiindPython(PyObject *module) {
    PyModule_AddFunctions(module, WinMiindPython_functions);

    PyModule_AddStringConstant(module, "__author__", "Hugh");
    PyModule_AddStringConstant(module, "__version__", "1.0.0");
    PyModule_AddIntConstant(module, "year", 2020);

    return 0; /* success */
}

/*
 * Documentation for WinMiindPython.
 */
PyDoc_STRVAR(WinMiindPython_doc, "The WinMiindPython module");


static PyModuleDef_Slot WinMiindPython_slots[] = {
    { Py_mod_exec, exec_WinMiindPython },
    { 0, NULL }
};

static PyModuleDef WinMiindPython_def = {
    PyModuleDef_HEAD_INIT,
    "WinMiindPython",
    WinMiindPython_doc,
    0,              /* m_size */
    NULL,           /* m_methods */
    WinMiindPython_slots,
    NULL,           /* m_traverse */
    NULL,           /* m_clear */
    NULL,           /* m_free */
};

PyMODINIT_FUNC PyInit_WinMiindPython() {
    return PyModuleDef_Init(&WinMiindPython_def);
}
