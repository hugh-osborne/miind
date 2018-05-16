# Once the Python Shared Library has been built in MIIND,
# copy this file to the results directory (where the .cpp and .so files were
# created).

import pylab
import numpy
import matplotlib.pyplot as plt
import imp

# Comment out MPI, comm and rank lines below if not using
# MPI
#######################
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
#######################

miind = imp.load_dynamic('libmiindlif', './libmiindlif.so')
number_of_nodes = 1
simulation_length = 3 #ms
miindmodel = miind.MiindModel(number_of_nodes, simulation_length)

miindmodel.init([])

# For MPI child processes, startSimulation runs the full simulation loop
# and so will not return until MPI process 0 has completed. At that point,
# we want to kill the child processes so that sim() is not called more than once.
if miindmodel.startSimulation() > 0 :
    quit()

constant_input = [5000]
activities = []
for i in range(int(simulation_length/0.001)): #0.001 is the time step defined in the xml
    activities.append(miindmodel.evolveSingleStep(constant_input)[0])

plt.figure()
plt.plot(activities)
plt.title("Firing Rate.")

plt.show()