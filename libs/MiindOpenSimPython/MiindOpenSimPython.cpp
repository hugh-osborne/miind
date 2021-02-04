#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "MiindOpenSim.h"

MiindOpenSim* model;

PyObject* miind_init(PyObject* self, PyObject* args)
{
    if (model) {
        delete model;
        model = NULL;
    }

    model = new MiindOpenSim();
    std::cout << "Miind object created.\n";
    model->buildSimulation();
    std::cout << "Built simulation.\n";

    Py_RETURN_NONE;
}

PyObject* miind_getTimeStep(PyObject* self, PyObject* args)
{
    return Py_BuildValue("d", model->getTimeStep());
}

PyObject* miind_getSimulationTime(PyObject* self, PyObject* args)
{
    return Py_BuildValue("d", model->getSimulationTime());
}

PyObject* miind_getCurrentSimulationTime(PyObject* self, PyObject* args)
{
    return Py_BuildValue("d", model->getCurrentSimulationTime());
}

PyObject* miind_startSimulation(PyObject* self, PyObject* args)
{
    model->beginSimulation();

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

    model->setInputs(activities);
    model->stepSimulation();

    Py_RETURN_NONE;
}

PyObject* miind_endSimulation(PyObject* self, PyObject* args)
{
    model->finishSimulation();

    Py_RETURN_NONE;
}

/*
 * List of functions to add to WinMiindPython in exec_WinMiindPython().
 */
static PyMethodDef MiindOpenSimPython_functions[] = {
    {"init",  miind_init, METH_VARARGS, "Init Miind Model."},
    {"getTimeStep",  miind_getTimeStep, METH_VARARGS, "Get time step."},
    {"getSimulationTime",  miind_getSimulationTime, METH_VARARGS, "Get sim time."},
    {"getCurrentSimulationTime", miind_getCurrentSimulationTime, METH_VARARGS, "Get current sim time."},
    {"startSimulation",  miind_startSimulation, METH_VARARGS, "Start simulation."},
    {"evolveSingleStep",  miind_evolveSingleStep, METH_VARARGS, "Evolve one time step."},
    {"endSimulation",  miind_endSimulation, METH_VARARGS, "Clean up."},
    { NULL, NULL, 0, NULL } /* marks end of array */
};

/*
 * Initialize WinMiindPython. May be called multiple times, so avoid
 * using static state. "oops." - Hugh
 */
int exec_MiindOpenSimPython(PyObject* module) {
    PyModule_AddFunctions(module, MiindOpenSimPython_functions);

    PyModule_AddStringConstant(module, "__author__", "Hugh");
    PyModule_AddStringConstant(module, "__version__", "1.0.0");
    PyModule_AddIntConstant(module, "year", 2020);

    return 0; /* success */
}

/*
 * Documentation for WinMiindPython.
 */
PyDoc_STRVAR(MiindOpenSimPython_doc, "The MiindOpenSimPython module");


static PyModuleDef_Slot MiindOpenSimPython_slots[] = {
    { Py_mod_exec, exec_MiindOpenSimPython },
    { 0, NULL }
};

static PyModuleDef MiindOpenSimPython_def = {
    PyModuleDef_HEAD_INIT,
    "MiindOpenSimPython",
    MiindOpenSimPython_doc,
    0,              /* m_size */
    NULL,           /* m_methods */
    MiindOpenSimPython_slots,
    NULL,           /* m_traverse */
    NULL,           /* m_clear */
    NULL,           /* m_free */
};

PyMODINIT_FUNC PyInit_MiindOpenSimPython() {
    return PyModuleDef_Init(&MiindOpenSimPython_def);
}
