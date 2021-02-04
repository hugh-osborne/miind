/* -------------------------------------------------------------------------- *
 *                 OpenSim:  Millard12EqMuscleWithAfferents.cpp               *
 * -------------------------------------------------------------------------- *
 */
 
//=============================================================================
// INCLUDES
//=============================================================================
#include "Millard12EqMuscleWithAfferents.h"
#include <iostream>  // remove later

//=============================================================================
// STATICS
//=============================================================================
using namespace std;
using namespace OpenSim;

///@cond  
const std::string Millard12EqMuscleWithAfferents::STATE_LPF_VELOCITY_NAME = "LPF_velocity";
const std::string Millard12EqMuscleWithAfferents::STATE_LPF_ACCELERATION_NAME = "LPF_acceleration";
// The below names are declared private in Millard2012EquilibriumMuscle so we have to duplicate them here. Boooo!
const string Millard12EqMuscleWithAfferents::AFF_STATE_ACTIVATION_NAME = "activation";
const string Millard12EqMuscleWithAfferents::AFF_STATE_FIBER_LENGTH_NAME = "fiber_length";
///@endcond  

//=============================================================================
// CONSTRUCTOR(S) AND DESTRUCTOR
//=============================================================================
/*
 * Default constructor.
 */
Millard12EqMuscleWithAfferents::Millard12EqMuscleWithAfferents()
{
	constructProperties();

	// Initialize the work variables to calculate the acceleration
	vel = 0.0;
	ts[4] = -0.004; ts[3] = -0.003; ts[2] = -0.002; ts[1] = -0.001; ts[0] = 0.0;
	C0 = 0.0; C1 = 0.0;
}

/*
 * Constructor.
 */
Millard12EqMuscleWithAfferents::Millard12EqMuscleWithAfferents(const std::string &name, 
					double maxIsometricForce, double optimalFiberLength, 
					double tendonSlackLength, double pennationAngle) : 
	Super(name, maxIsometricForce, optimalFiberLength, tendonSlackLength, 
          pennationAngle)
{
	constructProperties();
	// Initialize the work variables to calculate the acceleration
	vel = 0.0;
	ts[4] = -0.004; ts[3] = -0.003; ts[2] = -0.002; ts[1] = -0.001; ts[0] = 0.0;
	C0 = 0.0; C1 = 0.0;
}

Millard12EqMuscleWithAfferents::Millard12EqMuscleWithAfferents(const Millard2012EquilibriumMuscle& muscle)
: Super(muscle.getName(), muscle.getMaxIsometricForce(), muscle.getOptimalFiberLength(), muscle.getTendonSlackLength(), 
          muscle.getPennationAngleAtOptimalFiberLength())
{
	constructProperties();
	
	// Initialize the work variables to calculate the acceleration
	vel = 0.0;
	ts[4] = -0.004; ts[3] = -0.003; ts[2] = -0.002; ts[1] = -0.001; ts[0] = 0.0;
	C0 = 0.0; C1 = 0.0;
}

/*
 * Construct and initialize properties.
 * All properties are added to the property set. Once added, they can be
 * read in and written to files.
 */
void Millard12EqMuscleWithAfferents::constructProperties()
{
	setAuthors("Sergio Verduzco from code by Ajay Seth");
	constructProperty_lpf_tau(0.1); // LPF time constant
}

// Define new states and their derivatives in the underlying system
void Millard12EqMuscleWithAfferents::extendAddToSystem(SimTK::MultibodySystem& system) const
{
	// Allow Millard2012EquilibriumMuscle to add its states, cache, etc.
	// to the system
	Super::extendAddToSystem(system);
	
	// low-pass filtered state variables used to calculate derivatives 
	addStateVariable(STATE_LPF_VELOCITY_NAME); // fiber velocity
	addStateVariable(STATE_LPF_ACCELERATION_NAME); // fiber acceleration
}

void Millard12EqMuscleWithAfferents::extendInitStateFromProperties(SimTK::State& s) const
{
    Super::extendInitStateFromProperties(s);
	
	// I'll init the state, but not from properties
	setLPFvelocity(s, 0.0);
	setLPFacceleration(s, 0.0);
										
	// Initialize the work variables to calculate the acceleration
	vel = 0.0;
	ts[4] = -0.004; ts[3] = -0.003; ts[2] = -0.002; ts[1] = -0.001; ts[0] = 0.0;
	C0 = 0.0; C1 = 0.0;
}

void Millard12EqMuscleWithAfferents::extendSetPropertiesFromState(const SimTK::State& s)
{
    Super::extendSetPropertiesFromState(s);
}

void Millard12EqMuscleWithAfferents::extendConnectToModel(Model& aModel)
{
	Super::extendConnectToModel(aModel);
}

//--------------------------------------------------------------------------
// GET & SET Properties
//--------------------------------------------------------------------------
void Millard12EqMuscleWithAfferents::setLPFtau(double aLPFtau) {
	set_lpf_tau(aLPFtau);
}

//--------------------------------------------------------------------------
// GET & SET States and their derivatives
//--------------------------------------------------------------------------

double Millard12EqMuscleWithAfferents::getLPFvelocity(const SimTK::State& s) const {
	return getStateVariableValue(s, STATE_LPF_VELOCITY_NAME);
}	
void Millard12EqMuscleWithAfferents::setLPFvelocity(SimTK::State& s, double Velocity) const {
	setStateVariableValue(s, STATE_LPF_VELOCITY_NAME, Velocity);
}	
double Millard12EqMuscleWithAfferents::getLPFacceleration(const SimTK::State& s) const {
	return getStateVariableValue(s, STATE_LPF_ACCELERATION_NAME);
}
void Millard12EqMuscleWithAfferents::setLPFacceleration(SimTK::State& s, double Acceleration) const {
	setStateVariableValue(s, STATE_LPF_ACCELERATION_NAME, Acceleration);
}

//=============================================================================
// COMPUTATION
//=============================================================================
void Millard12EqMuscleWithAfferents::
computeInitialFiberEquilibrium(SimTK::State& s) const
{
	// First let the muscle find an equilibrium state
	Super::computeInitialFiberEquilibrium(s);
	
	setLPFvelocity(s, 0.0);
	// a simplifying assumption is a steady state
	setLPFacceleration(s, 0.0); 
	
	// update the work vectors assuming no acceleration
	vel = 0.0;
	ts[4] = -0.004; ts[3] = -0.003; ts[2] = -0.002; ts[1] = -0.001; ts[0] = 0.0;

	GTO.initFromMuscle(s);
}

void Millard12EqMuscleWithAfferents::computeStateVariableDerivatives(const SimTK::State& s) const
{	
	// vector of the derivatives to be returned
	SimTK::Vector derivs(getNumStateVariables(), 0.0);
	int nd = derivs.size();

	SimTK_ASSERT1(nd == 4, "Millard12EqMuscleWithAfferents: Expected 4 state variables"
        " but encountered  %f.", nd);

// This is the parent's computeStateVariableDerivatives
/*--------------------------------------------------------------------
	int idx = 0;

    if (!isDisabled(s)) {
        // Activation is the first state (if it is a state at all)
        if(!get_ignore_activation_dynamics() &&
           idx+1 <= getNumStateVariables()) {
               derivs[idx] = getActivationDerivative(s);
               idx++;
        }

        // Fiber length is the next state (if it is a state at all)
        if(!get_ignore_tendon_compliance() && idx+1 <= getNumStateVariables()) {
            derivs[idx] = getFiberVelocity(s);
        }
    }
--------------------------------------------------------------------*/
// This is a "carefree" version of that:
	derivs[0] = getActivationDerivative(s);
	derivs[1] = getFiberVelocity(s);
	
	// next state is the LPF velocity
	derivs[2] = (getFiberVelocity(s) - getLPFvelocity(s)) / 0.01;
	 
	// the LPF acceleration
	derivs[3] = (approxFiberAcceleration(s) - getLPFacceleration(s)) / 0.001;

	setStateVariableDerivativeValue(s, AFF_STATE_ACTIVATION_NAME, derivs[0]);
	setStateVariableDerivativeValue(s, AFF_STATE_FIBER_LENGTH_NAME, derivs[1]);
	setStateVariableDerivativeValue(s, STATE_LPF_VELOCITY_NAME, derivs[2]);
	setStateVariableDerivativeValue(s, STATE_LPF_ACCELERATION_NAME, derivs[3]);
}

//--------------------------------------------------------------------------
// Approximate the muscle fiber acceleration
//--------------------------------------------------------------------------
// HO : Originally approxFiberAcceleration used a more sophisticated way to estimate the
// acceleration but it doesn't play nice with a static time step which we need
// for syncing with the neural simulation so we just use the less accurate, naive, but
// faster solution just using the previous value.
//double Millard12EqMuscleWithAfferents::
//       approxFiberAcceleration(const SimTK::State& s) const
//{
//	double curr_time = s.getTime();
//	double curr_vel = getLPFvelocity(s);
//	double v;
//
//	if (curr_time - ts(0) > 0.0) {
//		v = (curr_vel - vel(0)) / (curr_time - ts(0));
//		ts(0) = curr_time;
//		vel(0) = curr_vel;
//	}
//	else {
//		v = 0.0;
//	}
//
//	return v;
//}

double Millard12EqMuscleWithAfferents::
approxFiberAcceleration(const SimTK::State& s) const
{
	// three point formula
	double curr_time = s.getTime();
	double curr_vel = getLPFvelocity(s);
	double v;

	if (curr_time - ts(0) > 0.0) {
		if (ts(3) >= 0.0) {
			v = ((-2 * vel(2)) + (9 * vel(1)) - (18 * vel(0)) + (11 * curr_vel)) / (6 * (curr_time - ts(0)));
			}
		else {
			v = 0.0;
		}

		ts(4) = ts(3); ts(3) = ts(2); ts(2) = ts(1); ts(1) = ts(0);
		ts(0) = curr_time;

		vel(4) = vel(3); vel(3) = vel(2); vel(2) = vel(1); vel(1) = vel(0);
		vel(0) = curr_vel;
	}
	else {
		v = C0(1,0);
	}
	C0(1, 0) = v;
	return v;
}