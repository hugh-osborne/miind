import os
import sys
import include
import algorithms
import nodes
import connections
import reporting
import simulation
import variables
import xml.etree.ElementTree as ET
import argparse
import directories3

XML_EXTENSION = '.xml'

# Nothing too fancy for the weight type

WEIGHTTYPES = ['double', 'DelayedConnection', 'CustomConnectionParameters']

def generate_preamble(outfile):
    outfile.write('//Machine-generated by miind.py. Edit at your own risk.\n\n')
    for inc in include.lib_includes:
        outfile.write(inc +'\n')

    return

def define_network_type(type):
    if type ==  "DelayedConnection":
        s = "MPILib::" + type
    elif type == "CustomConnectionParameters":
        s = "MPILib::" + type
    else:
        s = "double"
    return 'typedef MPILib::MPINetwork<' + s + ', MPILib::utilities::CircularDistribution> Network;\n'

def abstract_type(type):
    if type ==  "DelayedConnection":
        s = "MPILib::" + type
    elif type == "CustomConnectionParameters":
        s = "MPILib::" + type
    else:
        s = "double"
    return 'MPILib::MiindTvbModelAbstract<' + s + ', MPILib::utilities::CircularDistribution>'

def define_abstract_type(type):
    if type ==  "DelayedConnection":
        s = "MPILib::" + type
    elif type == "CustomConnectionParameters":
        s = "MPILib::" + type
    else:
        s = "double"
    return 'define_python_MiindTvbModelAbstract<' + s + ', MPILib::utilities::CircularDistribution>();\n'


def generate_closing(outfile, steps, t_step, weighttype, tree, prog_name):

    outfile.write('\t}\n')
    outfile.write('\t\n')
    outfile.write('\treport_handler = new MPILib::report::handler::InactiveReportHandler();\n')
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

    outfile.write('\n\nprivate:\n\n')
    outfile.write('\tunsigned long _count;\n')
    outfile.write('\tstd::vector<MPILib::NodeId> _display_nodes;\n');
    outfile.write('\tstd::vector<MPILib::NodeId> _rate_nodes;\n');
    outfile.write('\tstd::vector<MPILib::Time> _rate_node_intervals;\n');
    outfile.write('\tstd::vector<MPILib::NodeId> _density_nodes;\n');
    outfile.write('\tstd::vector<MPILib::Time> _density_node_start_times;\n');
    outfile.write('\tstd::vector<MPILib::Time> _density_node_end_times;\n');
    outfile.write('\tstd::vector<MPILib::Time> _density_node_intervals;\n');
    outfile.write('\n')
    for t in  algorithms.RATEFUNCTIONS:
        outfile.write(t)

    outfile.write('\n\n')

    outfile.write('};\n\n')

    type = weighttype.text
    variable_list = tree.findall('Variable')
    outfile.write('BOOST_PYTHON_MODULE(lib' + prog_name + ')\n')
    outfile.write('{\n')
    outfile.write('\tusing namespace boost::python;\n')
    outfile.write('\t' + define_abstract_type(type))
    outfile.write('\tclass_<MiindModel, bases<' + abstract_type(type) + '>>("MiindModel", init<int,long>())\n')

    if len(variable_list) > 0:
        var_types = variables.parse_variable_types(variable_list)
        outfile.write('\t.def(init<int,long' + var_types + '>())\n')
    outfile.write('\t.def("init", &MiindModel::init)\n')
    outfile.write('\t.def("init", &MiindModel::init)\n')
    outfile.write('\t.def("startSimulation", &MiindModel::startSimulation)\n')
    outfile.write('\t.def("endSimulation", &MiindModel::endSimulation)\n')
    outfile.write('\t.def("evolveSingleStep", &MiindModel::evolveSingleStep);\n')
    outfile.write('}\n')

def parse_xml(infile, outfile):
    tree=ET.fromstring(infile.read())
    m=tree.find('WeightType')
    s = m.text
    return define_network_type(s), tree

def constructor_override(outfile,tree,typ):

    variable_list = tree.findall('Variable')
    variables.parse_variables(variable_list,outfile)

    outfile.write('\tMiindModel(int num_nodes, long simulation_length ):\n')
    outfile.write('\t\tMiindTvbModelAbstract(num_nodes, simulation_length),_count(0){\n')
    outfile.write('#ifdef ENABLE_MPI\n')
    outfile.write('\t// initialise the mpi environment this cannot be forwarded to a class\n')
    outfile.write('\tboost::mpi::environment env();\n')
    outfile.write('#endif\n')
    outfile.write('}\n\n')

    outfile.write('\tMiindModel(long simulation_length ):\n')
    outfile.write('\t\tMiindTvbModelAbstract(1, simulation_length),_count(0){\n')
    outfile.write('#ifdef ENABLE_MPI\n')
    outfile.write('\t// initialise the mpi environment this cannot be forwarded to a class\n')
    outfile.write('\tboost::mpi::environment env();\n')
    outfile.write('#endif\n')
    outfile.write('}\n\n')

    if len(variable_list) > 0:
        outfile.write('\tMiindModel(int num_nodes, long simulation_length \n')
        variables.parse_variables_as_parameters(variable_list,outfile)
        outfile.write('):\n')
        outfile.write('\t\tMiindTvbModelAbstract(num_nodes, simulation_length),_count(0)\n')
        variables.parse_variables_as_constructor_defaults(variable_list, outfile)
        outfile.write('{\n')
        outfile.write('#ifdef ENABLE_MPI\n')
        outfile.write('\t// initialise the mpi environment this cannot be forwarded to a class\n')
        outfile.write('\tboost::mpi::environment env();\n')
        outfile.write('#endif\n')
        outfile.write('}\n\n')

    if len(variable_list) > 0:
        outfile.write('\tMiindModel(long simulation_length \n')
        variables.parse_variables_as_parameters(variable_list,outfile)
        outfile.write('):\n')
        outfile.write('\t\tMiindTvbModelAbstract(1, simulation_length),_count(0)\n')
        variables.parse_variables_as_constructor_defaults(variable_list, outfile)
        outfile.write('{\n')
        outfile.write('#ifdef ENABLE_MPI\n')
        outfile.write('\t// initialise the mpi environment this cannot be forwarded to a class\n')
        outfile.write('\tboost::mpi::environment env();\n')
        outfile.write('#endif\n')
        outfile.write('}\n\n')

def function_overrides(outfile,tree,typ):
    outfile.write('\n\tvoid endSimulation(){\n')
    outfile.write('\t\t'+ abstract_type(typ) +'::endSimulation();\n')
    outfile.write('\t}\n')

    node_list = tree.findall('Nodes/Node')
    nodemap = node_name_to_node_id(node_list)

    outfile.write('\tint startSimulation(){\n')
    outfile.write(reporting.define_display_nodes(tree,nodemap))
    outfile.write(reporting.define_rate_nodes(tree,nodemap))
    outfile.write(reporting.define_density_nodes(tree,nodemap))
    outfile.write('\n')
    outfile.write('\t\t_rate_nodes = rate_nodes;\n')
    outfile.write('\t\t_rate_node_intervals = rate_node_intervals;\n')
    outfile.write('\t\t_display_nodes = display_nodes;\n')
    outfile.write('\t\t_density_nodes = density_nodes;\n')
    outfile.write('\t\t_density_node_start_times = density_node_start_times;\n')
    outfile.write('\t\t_density_node_end_times = density_node_end_times;\n')
    outfile.write('\t\t_density_node_intervals = density_node_intervals;\n')
    outfile.write('\t\t\n')
    outfile.write('\t\tif (_display_nodes.size() > 0)\n')
    outfile.write('\t\t\tTwoDLib::Display::getInstance()->animate(true, display_nodes,_time_step);\n')
    outfile.write('\t\treturn '+ abstract_type(typ) +'::startSimulation();\n')
    outfile.write('\t}\n')

    t_step = tree.find('SimulationRunParameter/t_step')
    outfile.write('\t\tboost::python::list evolveSingleStep(boost::python::list c){\n')
    outfile.write('\t\tnetwork.reportNodeActivities(_rate_nodes, _rate_node_intervals, (_count * ' + t_step.text + '));\n')
    outfile.write('\t\tif (_display_nodes.size() > 0)\n')
    outfile.write('\t\t\tTwoDLib::Display::getInstance()->updateDisplay(_count);\n')
    outfile.write('\t\tTwoDLib::GridReport<'+typ+'>::getInstance()->reportDensity(_density_nodes,_density_node_start_times,_density_node_end_times,_density_node_intervals,(_count * _time_step));\n')
    outfile.write('\t\t_count++;\n')
    outfile.write('\t\treturn '+ abstract_type(typ) +'::evolveSingleStep(c);\n')
    outfile.write('\t}\n\n')

def generate_opening(outfile, tree, typ):
    outfile.write('class MiindModel : public ' + abstract_type(typ) + ' {\n')
    outfile.write('public:\n\n')
    constructor_override(outfile, tree,typ)
    function_overrides(outfile,tree,typ)
    outfile.write('\n')
    outfile.write('\tvoid init(boost::python::list params)\n')
    outfile.write('\t{\n')
    t_step = tree.find('SimulationRunParameter/t_step')
    outfile.write('\t\t_time_step = ' + t_step.text + ';\n')
    outfile.write('\t\tfor(int i=0; i<_num_nodes; i++) {\n')

def matrix_transform_name(fn):
    '''Identifies matrix transform files mentioned in an XML file. For example used in placing the right model file in
    the same directory as an XML file.'''
    infile = open(fn)
    tree=ET.fromstring(infile.read())
    ma = tree.findall('Algorithms/Algorithm')

    tmatnames = []
    for a in ma:
        if a.attrib['type'] in ['GridAlgorithm','GridAlgorithmGroup','GridJumpAlgorithm','GridSomaDendriteAlgorithm']:
            tmatnames.append(a.attrib['transformfile'])
    return tmatnames

def model_name(fn):
    '''Identifies model files mentioned in an XML file. For example used in placing the right model file in
    the same directory as an XML file.'''
    infile = open(fn)
    tree=ET.fromstring(infile.read())
    ma = tree.findall('Algorithms/Algorithm')

    modelnames = []
    for a in ma:
        if a.attrib['type'] in ['GridAlgorithm','GridAlgorithmGroup','GridJumpAlgorithm','GridSomaDendriteAlgorithm','MeshAlgorithm','MeshAlgorithmGroup']:
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

def node_name_to_node_id(nodes):
     '''Create a map from name to NodeId from node elements. Return this map.'''
     d ={}
     for i,node in enumerate(nodes):
          d[node.attrib['name']] = i
     return d

def generate_outputfile(infile, outfile, prog_name):
    generate_preamble(outfile)
    nettype, tree = parse_xml(infile,outfile)
    outfile.write(nettype)
    algies = tree.findall('Algorithms')
    if len(algies) != 1:
        raise ValueError

    alg_list = algies[0].findall('Algorithm')
    weighttype = tree.find('WeightType')
    generate_opening(outfile, tree, weighttype.text)
    outfile.write('\t// generating algorithms\n')
    algorithms.parse_algorithms(alg_list,weighttype,outfile,for_lib=True)
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

    t_begin = tree.find('SimulationRunParameter/t_begin')
    t_end   = tree.find('SimulationRunParameter/t_end')
    t_step = tree.find('SimulationRunParameter/t_step')

    generate_closing(outfile, '(' + t_end.text + ' - ' + t_begin.text + ') / ' + t_step.text , t_step.text, weighttype, tree, prog_name)

    algorithms.reset_algorithms()
    nodes.reset_nodes()
