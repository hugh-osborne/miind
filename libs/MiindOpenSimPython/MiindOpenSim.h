#ifndef INCLUDE_GUARD_MIINIDOPENSIM
#define INCLUDE_GUARD_MIINIDOPENSIM

#include "Millard12EqMuscleWithAfferents.h"
#include "NeuralController.h"

class MiindOpenSim {
public:
	MiindOpenSim();
	~MiindOpenSim();

	void buildSimulation();
	void beginSimulation();
	void stepSimulation();
	void finishSimulation();
	void setInputs(std::vector<double> inputs) { neural_controller->setInputOverrides(inputs); }

	double getSimulationTime() { return neural_controller->getSimulationTime(); }
	double getTimeStep() { return neural_controller->getSimulationTimestep(); }
	double getCurrentSimulationTime() { return neural_controller->getCurrentTime(); }

private:

	OpenSim::Model *model;
	OpenSim::TableReporter* table_reporter;
	OpenSim::ConsoleReporter* reporter;
	SimTK::RungeKutta3Integrator *integrator;
	SimTK::TimeStepper *time_stepper;
	OpenSim::NeuralController* neural_controller;
	void setElbowExtensionPosture();
};

#endif