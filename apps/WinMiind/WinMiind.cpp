#include <boost/timer/timer.hpp>
#include <MPILib/include/MPINetworkCode.hpp>
#include <MPILib/include/RateAlgorithmCode.hpp>
#include <MPILib/include/SimulationRunParameter.hpp>
#include <MPILib/include/WilsonCowanAlgorithm.hpp>
#include <MPILib/include/PersistantAlgorithm.hpp>
#include <MPILib/include/DelayAlgorithmCode.hpp>
#include <MPILib/include/RateFunctorCode.hpp>

#include "SimulationParser.h"

void main() {

	SimulationParser<MPILib::CustomConnectionParameters> sim_parser(std::string("lif.xml"));
	sim_parser.init();
	sim_parser.startSimulation();

	std::vector<double> inputs;
	inputs.push_back(0);

	while (!sim_parser.simulationComplete()) {
		sim_parser.evolveSingleStep(inputs);
	}

	std::cout << "complete.\n";
}