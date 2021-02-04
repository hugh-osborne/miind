#ifndef OPENSIM_MILEUSNIC_06_SPINDLE_H_
#define OPENSIM_MILEUSNIC_06_SPINDLE_H_
/*
* The Mileusnic06Spindle class implements the spindle afferent model from: 
* Mileusnic et al. "Mathematical Models of Proprioceptors. I. Control and
* Transduction in the Muscle Spinle" J Neurophysiol 96: 1772-1788, 2006.
* Objects of this class are meant to be added as subcomponents of muscle
* objects. Those muscles receive the inputs corresponding to static and dynamic 
* gamma motoneuron activity, in addition to their excitation. The fusimotor
* inputs have units of pulses per second.
*
* @author Sergio Verduzco Flores
*/

#include <OpenSim/Actuators/osimActuatorsDLL.h>
#include "OpenSim/Simulation/Model/ModelComponent.h"
#include "OpenSim/Simulation/Model/Model.h"
#include <OpenSim/Simulation/Model/Muscle.h>

namespace OpenSim {

class Mileusnic06Spindle : public ModelComponent {
OpenSim_DECLARE_CONCRETE_OBJECT(Mileusnic06Spindle, ModelComponent);
friend class Millard12EqMuscleWithAfferents;

public:
//==============================================================================
// PROPERTIES
//==============================================================================

	OpenSim_DECLARE_PROPERTY(default_activation, double,
		"default value for both static and dynamic activation");

	OpenSim_DECLARE_OUTPUT(primary_Ia, double,
            getIaOutput, SimTK::Stage::Dynamics);

	OpenSim_DECLARE_OUTPUT(secondary_II, double,
            getIIOutput, SimTK::Stage::Dynamics);

//=============================================================================
// METHODS
//=============================================================================

	Mileusnic06Spindle();  // the constructor 

	//-------------------------------------------------------------------------
	// GET & SET state variables and their derivatives
	//-------------------------------------------------------------------------
	
	double getDynamicActivation(const SimTK::State& s) const;
	void setDynamicActivation(SimTK::State& s, double Activation) const;
	double getStaticActivation(const SimTK::State& s) const;
	void setStaticActivation(SimTK::State& s, double Activation) const;
	
	double getTensionBag1(const SimTK::State& s) const;
	void setTensionBag1(SimTK::State& s, double Tension) const;
	double getTensionBag1Deriv(const SimTK::State& s) const;
	void setTensionBag1Deriv(SimTK::State& s, double TensionDeriv) const;
	
	double getTensionBag2(const SimTK::State& s) const;
	void setTensionBag2(SimTK::State& s, double Tension) const;
	double getTensionBag2Deriv(const SimTK::State& s) const;
	void setTensionBag2Deriv(SimTK::State& s, double TensionDeriv) const;
	
	double getTensionChain(const SimTK::State& s) const;
	void setTensionChain(SimTK::State& s, double Tension) const;
	double getTensionChainDeriv(const SimTK::State& s) const;
	void setTensionChainDeriv(SimTK::State& s, double TensionDeriv) const;
	
	double getIaOutput(const SimTK::State& s) const;
	double& updIaOutput(const SimTK::State& s) const;
	double getIIOutput(const SimTK::State& s) const;
	double& updIIOutput(const SimTK::State& s) const;
	
	//-------------------------------------------------------------------------
	// GET & SET Properties
	//-------------------------------------------------------------------------
	
	/** default static and dynamic activation  */
	double getDefaultActivation() const { return get_default_activation(); }
	void setDefaultActivation(double aDefaultActivation);
	
	//-------------------------------------------------------------------------
	// Compute initial conditions
	//-------------------------------------------------------------------------
	/** This method is called from the owner muscle's
	 *  computeInitialFiberEquilibrium method */
	void computeInitialSpindleEquilibrium(SimTK::State& s) const;
	
protected:

	static const std::string STATE_dynamic_activation_NAME;
	static const std::string STATE_static_activation_NAME;
	static const std::string STATE_tension_bag1_NAME;
	static const std::string STATE_tension_bag2_NAME;
	static const std::string STATE_tension_chain_NAME;
	static const std::string STATE_tension_bag1_deriv_NAME;
	static const std::string STATE_tension_bag2_deriv_NAME;
	static const std::string STATE_tension_chain_deriv_NAME;
	static const std::string CACHE_primaryIa_NAME;
	static const std::string CACHE_secondaryII_NAME;
	
	//--------------------------------------------------------------------------
	// MODEL COMPONENT INTERFACE
	//--------------------------------------------------------------------------
    /// Currently you need to call all these from the owner muscle.
	
	/** Add any Symbody elements(state variables, cache variables, etc.) 
	 * to the System. To be called after connectToModel calls. 
	 * Due to reasons not entirely clear, this method has to be invoked
	 * from the addToSystem method of the owner muscle. At the end. */
	void extendAddToSystem(SimTK::MultibodySystem& system) const OVERRIDE_11;
	
	/** initialize state variables from properties */
	void extendInitStateFromProperties(SimTK::State& s) const OVERRIDE_11;
	
	/** set current values of state into the properties */
	void extendSetPropertiesFromState(const SimTK::State& s) OVERRIDE_11;
	
	/** Extension of the parent class method. 
	 * This method must be called from the connectToModel method of
	 * the owner muscle, before the Super::connectToModel statement. */
	void extendConnectToModel(Model& aModel) OVERRIDE_11;
	
	//-------------------------------------------------------------------------
	// PARAMETERS
	//-------------------------------------------------------------------------	
		
	/** Not planning to alter the parameters from the Mileusnic06 paper, so to
	 *  to keep the code more manageable I didn't make the 52 parameters 
	 *  serializable. Instead, I put the parameters inside structs (one for each
	 *  type of fiber), and wrote struct constructors that filled all the 
	 *  default values. */
	struct bag1Params {
		double freq;	// pulses per second
		double p;		// dimensionless
		double tau;		// seconds
		double beta_0;	// force units / (L_0 / s) --> N s / m 
		double beta_1;	// force units / (L_0 / s) --> N s / m
		double Gamma_1;	// force units
		double C_L;		// dimensionless
		double C_S;		// dimensionless
		double R;		// length normalized by optimal fascicle length 
		double a; 		// dimensionless
		double K_SR;	// (force units) / L_0 
		double K_PR;	// (force units) / L_0
		double M;		// (force units) / (L_0/s^2) --> kg
		double L_0SR;	// L_0  (optimal fascicle length)
		double L_0PR;	// L_0 
		double L_NSR;	// L_0
		double G;		// (pulses per second)/(force units) 
		
		// this constructor initializes the parameters with the values from
		// TABLE 1 of Mileusnic 2006 I.
		bag1Params():
			tau(0.149),
			freq(60.0),
			p(2.0),
			beta_0(0.0605),
			beta_1(0.2592),
			Gamma_1(0.0289),
			C_L(1),
			C_S(0.42),
			R(0.46),
			a(0.3),
			K_SR(10.4649),
			K_PR(0.15),
			M(0.0002), 
			L_0SR(0.04),
			L_0PR(0.76),
			L_NSR(0.0423),
			G(20000) {}; // 20000 for the cat's spindle
	} bag1;
		
	struct bag2Params {
		double freq;	// pulses per second
		double p;		// dimensionless
		double tau;		// seconds
		double beta_0;	// (force units) / (L_0 / s) --> N s / (optimal length) 
		double beta_2;	// (force units) / (L_0 / s) --> N s / (optimal length)
		double Gamma_2;	// force units
		double C_L;		// dimensionless
		double C_S;		// dimensionless
		double R;		// length normalized by optimal fascicle length
		double a; 		// dimensionless
		double K_SR;	// (force units) / L_0 
		double K_PR;	// (force units) / L_0
		double M;		// (force units) / (L_0/s^2) --> kg
		double L_0SR;	// L_0  (optimal fascicle length)
		double L_0PR;	// L_0 
		double L_NSR;	// L_0
		double L_NPR;	// L_0
		double G_pri;	// pps/fu . When G contributes to primary afferent
		double G_sec;	// pps/fu . When G contributes to secondary afferent
		double X;		// dimensionless
		double L_sec;	// L_0
		
		// this constructor initializes the parameters with the values from
		// TABLE 1 of Mileusnic 2006 I.
		bag2Params():
			tau(0.205),
			freq(60.0),
			p(2.0),
			beta_0(0.0822),
			beta_2(-0.046),
			Gamma_2(0.0636),
			C_L(1),
			C_S(0.42),
			R(0.46),
			a(0.3),
			K_SR(10.4649),
			K_PR(0.15),
			M(0.0002), 
			L_0SR(0.04),
			L_0PR(0.76),
			L_NSR(0.0423),
			L_NPR(0.89),
			G_pri(10000),	// 10000 for the cat's spindle
			G_sec(7250),	// 7250 for the cat's spindle
			X(0.7),
			L_sec(0.04) {};
	} bag2;
	
	struct chainParams {
		double freq;	// pulses per second
		double p; 		// dimensionless
		double beta_0;	// (force units) / (L_0 / s) --> N s / (optimal length) 
		double beta_2;	// (force units) / (L_0 / s) --> N s / (optimal length)
		double Gamma_2;	// force units
		double C_L;		// dimensionless
		double C_S;		// dimensionless
		double R;		// length normalized by optimal fascicle length
		double a; 		// dimensionless
		double K_SR;	// (force units) / L_0 
		double K_PR;	// (force units) / L_0
		double M;		// (force units) / (L_0/s^2) --> kg
		double L_0SR;	// L_0  (optimal fascicle length)
		double L_0PR;	// L_0 
		double L_NSR;	// L_0
		double L_NPR;	// L_0
		double G_pri;	// pps/fu . When G contributes to primary afferent
		double G_sec;	// pps/fu . When G contributes to secondary afferent
		double X;		// dimensionless
		double L_sec;	// L_0
		
		// this constructor initializes the parameters with the values from
		// TABLE 1 of Mileusnic 2006 I.
		chainParams():
			freq(90.0),
			p(2.0),
			beta_0(0.0822),
			beta_2(-0.069),
			Gamma_2(0.0954),
			C_L(1),
			C_S(0.42),
			R(0.46),
			a(0.3),
			K_SR(10.4649),
			K_PR(0.15),
			M(0.0002), 
			L_0SR(0.04),
			L_0PR(0.76),
			L_NSR(0.0423),
			L_NPR(0.89),
			G_pri(10000),	// 10000 for the cat's spindle
			G_sec(7250),	// 7250 for the cat's spindle
			X(0.7),
			L_sec(0.04) {};
	} chain; 
	
	double S; // Amount of partial occlusion. S = 0.156 (Pg. 1776)
	
	//-------------------------------------------------------------------------
	// COMPUTATIONS
	//-------------------------------------------------------------------------

	void computeStateVariableDerivatives(const SimTK::State& s) const OVERRIDE_11;
		
private:

	void constructProperties();
	
	/** This pointer leads to the the owner muscle. 
	 *  It can't be used to call non-const functions.
	 *  Initialized in Mileusnic06Spindle::connectToModel */
	const Muscle *musclePtr;
	
	/** This is a utility function to get the sign of its argument.
	 *  Returns 1 if positive, -1 if negative, 0 otherwise. */
	int sgn(const double val) const {
		return ((0 < val) - (0 > val));
	}
	
	/** This is a utility function to get the max of two numbers */
	double max(double val1, double val2) const {
		return (val1 > val2)? val1 : val2;
	}

	/** This is a utility function to get the min of two numbers */
	double min(double val1, double val2) const {
		return (val1 < val2)? val1 : val2;
	}
	
//=============================================================================
};	// END of class Mileusnic06Spindle	
	
//=============================================================================
} // end of namespace OpenSim

#endif