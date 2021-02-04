# Once the Python Shared Library has been built in MIIND,
# copy this file to the results directory (where the .cpp and .so files were
# created).

import pylab
import numpy as np
import matplotlib.pyplot as plt
import csv

import sys
sys.path.insert(1, '../../x64/Release')
import MiindOpenSimPython as miind

try:
    miind.init()
except:
    print("init threw an exception but I don't know why and it doesn't apparently matter.")

timestep = miind.getTimeStep()
simulation_length = miind.getSimulationTime()
print('Timestep from XML : {}'.format(timestep))
print('Sim time from XML : {}'.format(simulation_length))

miind.startSimulation()

constant_input = [80000,180000,180000]
time = 0.0
while time < simulation_length:
    time = miind.getCurrentSimulationTime()
    if time > 1:
        if time > 3:
            constant_input = [80000,180000,180000]
        else:
            if time > 1.5:
                constant_input = [100000,180000,180000]
            else:
                constant_input = [80000+(((time-1.0)/0.5)*20000),180000,180000]
    miind.evolveSingleStep(constant_input)

miind.endSimulation()

rows = []

with open("output.txt") as outfile:
    filereader = csv.reader(outfile, delimiter=',')
    for line in filereader:
        ls = [float(l) for l in line[:-1]]
        rows = rows + [ls]

rows = np.array(rows)
rows = np.transpose(rows)
fig = plt.figure()

for i in range(7):
    norm = np.linalg.norm(rows[i])
    plt.plot(rows[i] / norm)

plt.legend(['Ia','II','Ib','Force','Length','Velocity','Acceleration'])
plt.show()
