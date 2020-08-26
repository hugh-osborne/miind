#include <boost/timer/timer.hpp>
#include <MPILib/include/MPINetworkCode.hpp>
#include <MPILib/include/RateAlgorithmCode.hpp>
#include <MPILib/include/SimulationRunParameter.hpp>
#include <MPILib/include/WilsonCowanAlgorithm.hpp>
#include <MPILib/include/PersistantAlgorithm.hpp>
#include <MPILib/include/DelayAlgorithmCode.hpp>
#include <MPILib/include/RateFunctorCode.hpp>

#include <WinMiindLib/SimulationParserCPU.h>

void main(int argc, char** argv) {
	if (argc == 1)
		std::cout << "WinMiind requires an XML simulation file.\n";
	if (argc > 2)
		std::cout << "WinMiind requires an XML simulation file only.\n";

	std::string current_exec_name = "coba.mmxml";
	SimulationParserCPU<MPILib::CustomConnectionParameters> sim_parser(current_exec_name);
	sim_parser.init();
	sim_parser.startSimulation();

	std::vector<double> inputs;
	inputs.push_back(0);

	while (!sim_parser.simulationComplete()) {
		sim_parser.evolveSingleStep(inputs);
	}

	std::cout << "complete.\n";
}