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
    for inc in include.lib_includes:
        outfile.write(inc +'\n')

    return

def define_network_type(type):
    if type ==  "DelayedConnection":
        s = "MPILib::" + type
    else:
        s = "double"
    return 'typedef MPILib::MPINetwork<' + s + ', MPILib::utilities::CircularDistribution> Network;\n'

def abstract_type(type):
    if type ==  "DelayedConnection":
        s = "MPILib::" + type
    else:
        s = "double"
    return 'MPILib::MiindTvbModelAbstract<' + s + ', MPILib::utilities::CircularDistribution>'

def define_abstract_type(type):
    if type ==  "DelayedConnection":
        s = "MPILib::" + type
    else:
        s = "double"
    return 'define_python_MiindTvbModelAbstract<' + s + ', MPILib::utilities::CircularDistribution>();\n'


def generate_closing(outfile, typ, tree):
    outfile.write('\t}\n')
    outfile.write('\t\n')
    name=tree.find('SimulationIO/SimulationName')
    name_str = name.text
    outfile.write('\tstd::string sim_name = \"' + name_str + '\";\n')
    outfile.write('#ifdef ENABLE_ROOT_REPORTER\n')
    outfile.write('\treport_handler = new MPILib::report::handler::RootReportHandler(sim_name,true);\n')
    outfile.write('#else\n')
    outfile.write('\treport_handler = new MPILib::report::handler::InactiveReportHandler();\n')
    outfile.write('#endif\n')
    outfile.write('\t\n')
    outfile.write('\tSimulationRunParameter par_run( *report_handler,(_simulation_length/_time_step)+1,0,\n')
    name=tree.find('SimulationRunParameter/name_log')
    log_name = name.text
    report_t_step = tree.find('SimulationRunParameter/t_report')
    report_state_t_step = tree.find('SimulationRunParameter/t_state_report')
    outfile.write('\t\t\t_simulation_length,' + report_t_step.text + ',_time_step,\"' + log_name + '\",'+ report_state_t_step.text + ');\n')
    outfile.write('\t\n')
    outfile.write('\tnetwork.configureSimulation(par_run);\n')
    outfile.write('\t}\n')
    outfile.write('};\n\n')

    for t in  algorithms.RATEFUNCTIONS:
        outfile.write(t)

    outfile.write('\n\n')
    variable_list = tree.findall('Variable')
    outfile.write('BOOST_PYTHON_MODULE(lib' + name_str + ')\n')
    outfile.write('{\n')
    outfile.write('\tusing namespace boost::python;\n')
    outfile.write('\t' + define_abstract_type(typ))
    outfile.write('\tclass_<MiindModel, bases<' + abstract_type(typ) + '>>("MiindModel", init<int,long>())\n')

    if len(variable_list) > 0:
        var_types = variables.parse_variable_types(variable_list)
        outfile.write('\t.def(init<int,long' + var_types + '>())\n')
    outfile.write('\t.def("init", &MiindModel::init);\n')
    outfile.write('}\n')

def parse_xml(infile, outfile):
    tree=ET.fromstring(infile.read())
    m=tree.find('WeightType')
    s = m.text
    return define_network_type(s), tree

def constructor_override(outfile,tree):

    variable_list = tree.findall('Variable')
    variables.parse_variables(variable_list,outfile)

    outfile.write('\tMiindModel(int num_nodes, long simulation_length ):\n')
    outfile.write('\t\tMiindTvbModelAbstract(num_nodes, simulation_length){}\n\n')

    if len(variable_list) > 0:
        outfile.write('\tMiindModel(int num_nodes, long simulation_length \n')
        variables.parse_variables_as_parameters(variable_list,outfile)
        outfile.write('):\n')
        outfile.write('\t\tMiindTvbModelAbstract(num_nodes, simulation_length)\n')
        variables.parse_variables_as_constructor_defaults(variable_list, outfile)
        outfile.write('{}\n\n')

def generate_opening(outfile, tree, typ):
    outfile.write('class MiindModel : public ' + abstract_type(typ) + ' {\n')
    outfile.write('public:\n\n')
    constructor_override(outfile, tree)
    outfile.write('\n')
    outfile.write('\tvoid init(boost::python::list params)\n')
    outfile.write('\t{\n')
    t_step = tree.find('SimulationRunParameter/t_step')
    outfile.write('\t\t_time_step = ' + t_step.text + ';\n')
    outfile.write('\t\tfor(int i=0; i<_num_nodes; i++) {\n')


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

    algies = tree.findall('Algorithms')
    if len(algies) != 1:
        raise ValueError

    alg_list = algies[0].findall('Algorithm')
    weighttype = tree.find('WeightType')
    generate_opening(outfile, tree, weighttype.text)
    outfile.write('\t// generating algorithms\n')
    algorithms.parse_algorithms(alg_list,weighttype,outfile)
    node_list = tree.findall('Nodes/Node')
    outfile.write('\t// generating nodes\n')
    nodes.parse_nodes(node_list,weighttype,outfile)
    outfile.write('\t// generating connections\n')
    connection_list = tree.findall('Connections/Connection')
    connections.parse_connections(connection_list,weighttype,outfile)
    connection_list = tree.findall('Connections/IncomingConnection')
    connections.parse_incoming_connections(connection_list,weighttype,outfile)
    connection_list = tree.findall('Connections/OutgoingConnection')
    connections.parse_outgoing_connections(connection_list,outfile)
    generate_closing(outfile, weighttype.text, tree)
    algorithms.reset_algorithms()
    nodes.reset_nodes()
