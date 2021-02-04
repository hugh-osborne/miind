#ifndef OPENSIM_GOLGI_TENDON_ORGAN_H_
#define OPENSIM_GOLGI_TENDON_ORGAN_H_

#include "OpenSim/Simulation/Model/ModelComponent.h"
#include "OpenSim/Simulation/Model/Model.h"
#include <OpenSim/Simulation/Model/Muscle.h>

namespace OpenSim {
	
/* 
* The Lin02GolgiTendonOrgan class implements a model for the aggregate
* population response of the Golgi tendon organs in a muscle.
* The model is taken from: Lin & Crago, "Neural and Mechanical Contributions
* to the Stretch Reflex: A Model Synthesis" Ann Biomed Eng 30:54-67.
*
* Objects of this class are meant to be contained in particular muscle
* objects, such as those of the Millard12EqMuscleWithAfferents class.
*
* @author Sergio Verduzco Flores
*/
class Lin02GolgiTendonOrgan : public ModelComponent {
OpenSim_DECLARE_CONCRETE_OBJECT(Lin02GolgiTendonOrgan, ModelComponent);
friend class Millard12EqMuscleWithAfferents;
friend class AfferentAnalysis;

public:
	OpenSim_DECLARE_PROPERTY(lpf_tau, double, 
        "time constant for the low-pass filters");

	OpenSim_DECLARE_OUTPUT(gto_out, double,
            getGTOout, SimTK::Stage::Dynamics);

	Lin02GolgiTendonOrgan();  // the constructor
	
	//-------------------------------------------------------------------------
	// GET & SET state variables and cache variables
	//-------------------------------------------------------------------------
	// low-pass filtered output of the nonlinearity
	double getX(const SimTK::State& s) const;
	void setX(SimTK::State& s, double X) const;
	// derivative of the low-pass filtered output of the nonlinearity
	double getXp(const SimTK::State& s) const;
	void setXp(SimTK::State& s, double Xp) const;
	// Intermediate variable for the transfer function filter
	double getY(const SimTK::State& s) const;
	void setY(SimTK::State& s, double Y) const;
	// output variable of the transfer function
	double getZ(const SimTK::State& s) const;
	void setZ(SimTK::State& s, double Z) const;
	// output of the Golgi tendon organ
	double getGTOout(const SimTK::State& s) const;
	double& updGTOout(const SimTK::State& s) const;
	
	//-------------------------------------------------------------------------
	// GET & SET Properties
	//-------------------------------------------------------------------------
	double getLPFtau() const { return get_lpf_tau(); }
	void setLPFtau(double aLPFtau);
	
protected:

	static const std::string STATE_LPF_OUTPUT_NAME;
	static const std::string STATE_LPF_DERIV_NAME;
	static const std::string STATE_FILTER_INTER_NAME;
	static const std::string STATE_FILTER_OUTPUT_NAME;

	static const std::string CACHE_GTO_OUT_NAME;

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
	
	/** Extension of the parent class method. 
	 * This method must be called from the connectToModel method of
	 * the owner muscle, before the Super::connectToModel statement. */
	void extendConnectToModel(Model& aModel) OVERRIDE_11;	
	
	//-------------------------------------------------------------------------
	// PARAMETERS
	//-------------------------------------------------------------------------	
	// Values used in Eq. 1 (initialized in the constuctor) 
	double Gg, Gf; // pulses/s and Newtons, respectively
	// Threshold for Eq. 3
	double thr;  // pulses / s
	
	/** Set an inital value for the output of the nonlinearity,
	  * and set the other variables to zero */
	void initFromMuscle(SimTK::State& s) const;
	
	//-------------------------------------------------------------------------
	// COMPUTATIONS
	//-------------------------------------------------------------------------
	void computeStateVariableDerivatives(const SimTK::State& s) const OVERRIDE_11;
	
	/** This method calculates the first and second derivatives of 
	 *  the log nonlinearity. 
	 *  It uses a version of the method in: Fornberg 1998 "Calculation of 
	 *  Weights in Finite Difference Formulas" SIAM Rev.40(3):685-691,
	 *  as well as simpler rules.
	 *  The stored data points to apply the method are in the vectors
	 *  ts, vel, F, and dF .  */
	SimTK::Vec<2> calculateDerivatives(const SimTK::State& s) const;
	
	
	
private:

	void constructProperties();
	
	/** The name of the muscle that owns the instance of this object.
	 *  Initialized in the muscle's connectToModel method, or in its constructor.
	 *  Checked in addToSystem. */
	std::string ownerMuscleName;
	
	/** This pointer leads to the the owner muscle. 
	 *  Initialized in Lin02GolgiTendonOrgan::connectToModel */
	const Muscle *musclePtr;
	
	/** These auxiliary vectors are used to approximate the 
	 *  derivatives of the log nonlinearity. 
	 *  'nl' and 'Dnl' are the stored values of the output of the
	 *  nonlinearity (Eq. 1) and its derivative, respectively. 
	 *  'ts' contains the time when they were stored.
	 *  Since computeStateVariableDerivatives is const, mutable is needed. */
	mutable SimTK::Vec<5> ts, nl, Dnl;
	mutable SimTK::Vec<2> previous_diffs;

}; // end of class Lin02GolgiTendonOrgan
	
} // end of namespace OpenSim

#endif