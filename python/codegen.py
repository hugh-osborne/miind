import os
import sys
import include
import algorithms
import nodes
import connections
import simulation
import variables
import xml.etree.ElementTree as ET
import argparse
import directories

XML_EXTENSION = '.xml'

# Nothing too fancy for the weight type

WEIGHTTYPES = ['double', 'DelayedConnection']

def generate_preamble(outfile):
    outfile.write('//Machine-generated by miind.py. Edit at your own risk.\n\n')
    for inc in include.includes:
        outfile.write(inc +'\n')

    return

def generate_closing(outfile):
    outfile.write('\tnetwork.configureSimulation(par_run);\n')
    outfile.write('\tnetwork.evolve();\n')
    outfile.write('\t} catch(std::exception& exc){\n')
    outfile.write('\t\tstd::cout << exc.what() << std::endl;\n')

    outfile.write('#ifdef ENABLE_MPI\n')
    outfile.write('\t//Abort the MPI environment in the correct way :\n')
    outfile.write('\tenv.abort(1);\n')
    outfile.write('#endif\n')
    outfile.write('\t}\n\n')

    outfile.write('\tMPILib::utilities::MPIProxy().barrier();\n')
    outfile.write('\tt.stop();\n')
    outfile.write('\tif (MPILib::utilities::MPIProxy().getRank() == 0) {\n')
    outfile.write('\n\t\tstd::cout << \"Overall time spend\\n\";\n')
    outfile.write('\t\tt.report();\n')
    outfile.write('\t}\n')

    outfile.write('\treturn 0;\n}\n')

    for t in  algorithms.RATEFUNCTIONS:
        outfile.write(t)

    return


def define_network_type(outfile, type):
    if type ==  "DelayedConnection":
        s = "MPILib::" + type
    else:
        s = "double"
    return 'typedef MPILib::MPINetwork<' + s + ', MPILib::utilities::CircularDistribution> Network;\n'

def parse_xml(infile, outfile):
    tree=ET.fromstring(infile.read())
    m=tree.find('WeightType')
    s = m.text
    return define_network_type(outfile,s), tree

def generate_opening(outfile):
    outfile.write('int main(int argc, char *argv[]){\n\tNetwork network;\n')
    outfile.write('\tboost::timer::auto_cpu_timer t;\n\n')
    outfile.write('#ifdef ENABLE_MPI\n')
    outfile.write('\t// initialise the mpi environment this cannot be forwarded to a class\n')
    outfile.write('\tboost::mpi::environment env(argc, argv);\n')
    outfile.write('#endif\n\n')
    outfile.write('\ttry {')


def model_name(fn):
    '''Identifies model files mentioned in an XML file. For example used in placing the right model file in
    the same directory as an XML file.'''
    infile = open(fn)
    tree=ET.fromstring(infile.read())
    ma = tree.findall('Algorithms/Algorithm')

    modelnames = []
    for a in ma:
        if a.attrib['type'] == 'MeshAlgorithm':
            modelnames.append(a.attrib['modelfile'])
    return modelnames

def matrix_names(fn):
    '''Find the file names of all MatrixFiles, mentioned in an XML file.'''
    infile = open(fn)
    tree=ET.fromstring(infile.read())
    ma = tree.findall('Algorithms/Algorithm/MatrixFile')
    matrixnames = []
    for a in ma:
        matrixnames.append(a.text)
    return matrixnames

def generate_outputfile(infile, outfile):
    generate_preamble(outfile)
    nettype, tree = parse_xml(infile,outfile)
    outfile.write(nettype)
    outfile.write('\t// defining variables\n') # whatever variables are use are global
    variable_list = tree.findall('Variable')
    variables.parse_variables(variable_list,outfile)
    algies = tree.findall('Algorithms')
    if len(algies) != 1:
        raise ValueError

    alg_list = algies[0].findall('Algorithm')
    weighttype = tree.find('WeightType')
    generate_opening(outfile)
    outfile.write('\t// generating algorithms\n')
    algorithms.parse_algorithms(alg_list,weighttype,outfile)
    node_list = tree.findall('Nodes/Node')
    outfile.write('\t// generating nodes\n')
    nodes.parse_nodes(node_list,weighttype,outfile)
    outfile.write('\t// generating connections\n')
    connection_list = tree.findall('Connections/Connection')
    connections.parse_connections(connection_list,weighttype,outfile)
    outfile.write('\t// generation simulation parameter\n')
    simhand = tree.find('SimulationIO')
    simulation.parse_simulation(simhand,outfile)
    simpar = tree.find('SimulationRunParameter')
    simulation.parse_parameter(simpar,outfile)
    generate_closing(outfile)
    algorithms.reset_algorithms()
    nodes.reset_nodes()
