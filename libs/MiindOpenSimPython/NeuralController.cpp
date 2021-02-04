#include "NeuralController.h"
#include <OpenSim/OpenSim.h>

namespace OpenSim {

	NeuralController::NeuralController(std::string xmlfile, const std::string& controlsFileName,
		int interpMethodType = 1)
		: Controller(),
		sim(new SimulationParserCPU<MPILib::CustomConnectionParameters>(xmlfile)) {
		constructProperties();
		sim->init();
	}

	NeuralController::NeuralController(std::string xmlfile)
		: Controller(),
		sim(new SimulationParserCPU<MPILib::CustomConnectionParameters>(xmlfile)) {
		constructProperties();
		sim->init();
	}

	NeuralController::~NeuralController() {}

	void NeuralController::setNull()
	{
		setAuthors("Hugh Osborne from code by Ajay Seth");
	}

	void NeuralController::constructProperties()
	{
		constructProperty_InputFunctions(FunctionSet());
	}

	void NeuralController::extendConnectToModel(Model& model)
	{
		Super::extendConnectToModel(model);
		// Probably would be better to call this in some init function which
		// I think there is but can't be bothered to check.
		sim->startSimulation();
	}

	void NeuralController::computeControls(const SimTK::State& s, SimTK::Vector& controls) const 
	{
		// Calculate the functions which define the inputs
		SimTK::Vector time(1, s.getTime());
		std::vector<double> inputs;
		for (int i = 0; i < get_InputFunctions().getSize(); i++)
			inputs.push_back(get_InputFunctions()[i].calcValue(time));

		if (input_overrides.size() > 0)
			inputs = input_overrides;

		// Run the simulation up to the current time so everthing's in sync
		std::vector<double> outputs;
		// Update MIIND and load the outputs into the functions
		while (s.getTime() > sim->getCurrentSimTime())
			outputs = sim->evolveSingleStep(inputs);
		
		if (outputs.size() <= 0)
			return;

		for (int i = 0; i < getActuatorSet().getSize(); i++) {
			SimTK::Vector actControls(3, 0.0);
			actControls[0] = outputs[sim->getIndexOfOutputNode(alpha_mapping[i])] / 50.0;
			actControls[1] = outputs[sim->getIndexOfOutputNode(beta_mapping[i])];
			actControls[2] = outputs[sim->getIndexOfOutputNode(gamma_mapping[i])];
			getActuatorSet()[i].addInControls(actControls, controls);
		}
	}

	void NeuralController::prescribeInputForControl(int input_index, Function* function) {
		if (input_index >= get_InputFunctions().getSize())
		{
			upd_InputFunctions().setSize(input_index + 1);
		}
		upd_InputFunctions().set(input_index, function);
	}

	void NeuralController::prescribeControlForActuator(int index, std::string alpha_node, std::string beta_node, std::string gamma_node) {

		if (index >= alpha_mapping.size())
		{
			alpha_mapping.resize(index + 1);
		}
		alpha_mapping[index] = alpha_node;

		if (index >= beta_mapping.size())
		{
			beta_mapping.resize(index + 1);
		}
		beta_mapping[index] = beta_node;

		if (index >= gamma_mapping.size())
		{
			gamma_mapping.resize(index + 1);
		}
		gamma_mapping[index] = gamma_node;

	}

	void NeuralController::prescribeControlForActuator(const std::string actName, std::string alpha_node, std::string beta_node, std::string gamma_node) {

		int index = getProperty_actuator_list().findIndex(actName);
		if (index < 0)
			throw Exception("NeuralController does not have " + actName + " in its list of actuators to control.");

		if (index >= alpha_mapping.size())
		{
			alpha_mapping.resize(index + 1);
		}
		alpha_mapping[index] = alpha_node;

		if (index >= beta_mapping.size())
		{
			beta_mapping.resize(index + 1);
		}
		beta_mapping[index] = beta_node;

		if (index >= gamma_mapping.size())
		{
			gamma_mapping.resize(index + 1);
		}
		gamma_mapping[index] = gamma_node;
	}
}

