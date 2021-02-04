#ifndef OPENSIM_MILLARD_EQ_MUSCLE_WITH_AFFERENTS_H_
#define OPENSIM_MILLARD_EQ_MUSCLE_WITH_AFFERENTS_H_

// INCLUDE
#include <OpenSim/OpenSim.h>
#include "Mileusnic06Spindle.h"
#include "Lin02GolgiTendonOrgan.h"

namespace OpenSim {

//=============================================================================
//=============================================================================
/**
 * This class extends a Millard2012EquilibriumMuscle by including including
 * inputs and outputs to represent spindle afferents and Golgi tendon organs. 
 * The spindle afferent model comes from: Mileusnic et al 2006 "Mathematical 
 * models of proprioceptors. I. Control and transduction in the muscle spindle"
 * J Neurophysiol 96:1772-1788.
 * The GTO model is the one used in: Lin, Crago 2002 "Neural and mechanical 
 * contributions to the stretch reflex: A model synthesis" 
 * Ann Biomed Eng 30:54-67.
 *
 * The implementation of this class was achieved by modifying the file
 * FatigableMuscle.h, which is part of the OpenSim API examples.
 *
 * @author Sergio Verduzco Flores (based on Millard2012EquilibriumMuscle)
 *
 * The Muscle base class specifies the interface that must be implemented 
 * by derived muscle classes. The Millard12EqMuscleWithAfferents derives from
 * Millard2012EquilibriumMuscle, which is a concrete implementation of the 
 * Muscle interface. The dynamics for spindle afferents are handled by
 * a member object of the Mileusnic06Spindle class. The dynamics of the
 * GTO are handled by an object of the Lin02GolgiTendonOrgan class.
 * 
 * The Mileusnic06Spindle object requires muscle fiber acceleration, which is 
 * not available from the Millard2012EquilibriumMuscle class, so this has to
 * be approximated numerically. To aid this, this muscle class provides a
 * low-pass filtered version of the fiber velocity, and the acceleration is 
 * provided as a state variable.
 *
 * @see Millard2012EquilibriumMuscle
 * @see Muscle
 * @see Mileusnic06Spindle
 * @see Lin02GolgiTendonOrgan  
 */
class Millard12EqMuscleWithAfferents : public Millard2012EquilibriumMuscle {
OpenSim_DECLARE_CONCRETE_OBJECT(Millard12EqMuscleWithAfferents, 
								Millard2012EquilibriumMuscle);
friend class Mileusnic06Spindle;
friend class Lin02GolgiTendonOrgan;

public:
//=============================================================================
// PROPERTIES
//=============================================================================
	OpenSim_DECLARE_PROPERTY(lpf_tau, double, 
        "time constant for all the low-pass filters");

	OpenSim_DECLARE_OUTPUT(lpf_velocity, double,
		getLPFvelocity, SimTK::Stage::Dynamics);

	OpenSim_DECLARE_OUTPUT(lpf_acceleration, double,
		getLPFacceleration, SimTK::Stage::Dynamics);

//=============================================================================
// THE AFFERENTS >><< >><< >><< >><< >><< >><< >><< >><< >><< >><< >><< >><<
//=============================================================================
	
	Mileusnic06Spindle spindle;
	Lin02GolgiTendonOrgan GTO;

private:

	/** These auxiliary vectors are used to approximate the 
	 *  derivatives of muscle variables. 
	 *  'vel' and 'ts' are the velocity and time values used in the method
	 *  to calculate the acceleration. 
	 *  Since computeStateVariableDerivatives is const, mutable is needed. */
	mutable SimTK::Vec<5> vel, ts;
	mutable SimTK::Mat33 C0;
	mutable SimTK::Mat44 C1;
		
public:
//=============================================================================
// METHODS
//=============================================================================
	/** The class constructors */
	Millard12EqMuscleWithAfferents();
	Millard12EqMuscleWithAfferents(const std::string &name, double maxIsometricForce, 
					double optimalFiberLength, double tendonSlackLength,
					double pennationAngle);
	Millard12EqMuscleWithAfferents(const Millard2012EquilibriumMuscle&);

    // employs the default destructor, copy constructor and copy assignment 
	// that are automatically supplied by the compiler if none are defined
	
	//-------------------------------------------------------------------------
	// GET & SET state variables and their derivatives
	//-------------------------------------------------------------------------
	double getLPFvelocity(const SimTK::State& s) const;
	void setLPFvelocity(SimTK::State& s, double Velocity) const;
	
	double getLPFacceleration(const SimTK::State& s) const;
	void setLPFacceleration(SimTK::State& s, double Acceleration) const;
	
	//-------------------------------------------------------------------------
	// GET & SET Properties
	//-------------------------------------------------------------------------
	double getLPFtau() const { return get_lpf_tau(); }
	void setLPFtau(double aLPFtau);
	
	//-------------------------------------------------------------------------
	// Other methods
	//-------------------------------------------------------------------------
	/** This function should override Actuator::numControls() */
	int numControls() const {return 3;};

	void setupAfferents() {
		addComponent(&spindle);
		addComponent(&GTO);
	}

	void setStageForOutputs(const Model& aModel, const SimTK::State& s) {
		if(!(getSpindle()->isCacheVariableValid(s, Mileusnic06Spindle::CACHE_primaryIa_NAME)) ||
		!(getSpindle()->isCacheVariableValid(s, Mileusnic06Spindle::CACHE_secondaryII_NAME)) ||
		!(getGTO()->isCacheVariableValid(s, Lin02GolgiTendonOrgan::CACHE_GTO_OUT_NAME)) )
		{
			aModel.getMultibodySystem().realize(s,SimTK::Stage::Acceleration);
		}
	}

	/** This function allows to peek at the Mileusnic06Spindle object.
	 *  Mostly for the sake of Analysis objects. */
	const Mileusnic06Spindle* getSpindle() { return &spindle; };
	
	/** This one allows to peek at the Lin02GolgiTendonOrgan object. */
	const Lin02GolgiTendonOrgan* getGTO() { return &GTO; };
	
protected:

	static const std::string STATE_LPF_VELOCITY_NAME;
	static const std::string STATE_LPF_ACCELERATION_NAME;
	static const std::string AFF_STATE_ACTIVATION_NAME;
	static const std::string AFF_STATE_FIBER_LENGTH_NAME;

	// Model Component Interface
	/** add new dynamical states to the multibody system corresponding
	    to this muscle */
	void extendAddToSystem(SimTK::MultibodySystem& system) const OVERRIDE_11;
	/** initialize muscle state variables from properties. For example, any 
	    properties that contain default state values */
	void extendInitStateFromProperties(SimTK::State& s)const OVERRIDE_11;
	/** use the current values in the state to update any properties such as 
	    default values for state variables */
	void extendSetPropertiesFromState(const SimTK::State& s) OVERRIDE_11;
	/** Extension of the parent class method.
	    Used here to declare the spindle as a subcomponent. */
	void extendConnectToModel(Model& aModel) OVERRIDE_11;

	//-------------------------------------------------------------------------
	// COMPUTATIONS
	//-------------------------------------------------------------------------
	/** This function finds an initial state for the tension, and
	 *  then it calls the homonimous method of the parent class. */
	void computeInitialFiberEquilibrium(SimTK::State& s) const OVERRIDE_11;
	
	/** Compute the derivatives for state variables added by this muscle */
	void computeStateVariableDerivatives(const SimTK::State& s) const OVERRIDE_11;
		
	/** The following function calculates the derivative of the fiber velocity
	 *  using the method in: Fornberg 1998 "Calculation of Weights in Finite
	 *  Difference Formulas" SIAM Rev.40(3):685-691.
	 *  The method is adapted to use with four points and to always evaluate at the
	 *  last point.  */
	double approxFiberAcceleration(const SimTK::State& s) const;
	double approxFiberVelocity(const SimTK::State& s) const;
	
private:

	/** construct the new properties and set their default values */
	void constructProperties();
	
public:
	/** This storage object records the afferent outputs */
	Storage *afferents;
	
//=============================================================================
};	// END of class Millard12EqMuscleWithAfferents
//=============================================================================
//=============================================================================

} // end of namespace OpenSim

#endif