
//=============================================================================
// INCLUDES
//=============================================================================
#include <iostream>  // for warning mesages
#include <cmath>  // for the pow() function
#include <stdlib.h> // for the abs function
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/ForceSet.h>
#include "Millard12EqMuscleWithAfferents.h"

//=============================================================================
// STATICS
//=============================================================================
using namespace std;
using namespace OpenSim;
using namespace SimTK;

///@cond  
const std::string Mileusnic06Spindle::STATE_dynamic_activation_NAME = "dynamic_activation";
const std::string Mileusnic06Spindle::STATE_static_activation_NAME = "static_activation";
const std::string Mileusnic06Spindle::STATE_tension_bag1_NAME = "tension_bag1";
const std::string Mileusnic06Spindle::STATE_tension_bag2_NAME = "tension_bag2";
const std::string Mileusnic06Spindle::STATE_tension_chain_NAME = "tension_chain";
const std::string Mileusnic06Spindle::STATE_tension_bag1_deriv_NAME = "tension_bag1_deriv";
const std::string Mileusnic06Spindle::STATE_tension_bag2_deriv_NAME = "tension_bag2_deriv";
const std::string Mileusnic06Spindle::STATE_tension_chain_deriv_NAME = "tension_chain_deriv";
const std::string Mileusnic06Spindle::CACHE_primaryIa_NAME = "primaryIa";
const std::string Mileusnic06Spindle::CACHE_secondaryII_NAME = "secondaryII";
///@endcond  

//=============================================================================
// CONSTRUCTOR(S) AND DESTRUCTOR
//=============================================================================

Mileusnic06Spindle::Mileusnic06Spindle() 
{
	constructProperties();

	// one parameter value that I didn't set anywhere else
	S = 0.156;
}

//=============================================================================
// MODEL COMPONENT INTERFACE
//=============================================================================

void Mileusnic06Spindle::extendAddToSystem(SimTK::MultibodySystem& system) const
{
	Super::extendAddToSystem(system);
	// adding state variables 
	addStateVariable(STATE_dynamic_activation_NAME);
	addStateVariable(STATE_static_activation_NAME);
	addStateVariable(STATE_tension_bag1_NAME);
	addStateVariable(STATE_tension_bag2_NAME);
	addStateVariable(STATE_tension_chain_NAME);
	addStateVariable(STATE_tension_bag1_deriv_NAME); // the first derivatives are also
	addStateVariable(STATE_tension_bag2_deriv_NAME); // state variables
	addStateVariable(STATE_tension_chain_deriv_NAME);

	// adding the output in cache variables
	addCacheVariable(CACHE_primaryIa_NAME, 0.0, SimTK::Stage::Dynamics);
	addCacheVariable(CACHE_secondaryII_NAME, 0.0, SimTK::Stage::Dynamics);
	
	// ensuring that the owner muscle is in the system	
	ForceSet& fSet = _model->updForceSet();

	auto ownerMuscleName = traversePathToComponent<Muscle>(getAbsolutePath().getParentPath())->getName();
	try {
		fSet.get(ownerMuscleName);
	}
	catch (OpenSim::Exception e) {
		std::cout << "WARNING - Mileusnic06Spindle::extendAddToSystem() could not find ";
		std::cout << "the muscle with name" << ownerMuscleName << '\n';
		std::cout << "Exception: " << e.getMessage() << '\n';
		return;
	}
	
	// ensuring that the owner muscle is of the right type
	std::string forceClassName = fSet.get(ownerMuscleName).getConcreteClassName();
	if( forceClassName != "Millard12EqMuscleWithAfferents" )
	{
		std::cout << "WARNING - In Mileusnic06Spindle::extendAddToSystem() \n";
		std::cout << "Mileusnic06Spindle is owned by a force that is not " ;
		std::cout << "of the Millard12EqMuscleWithAfferents class \n" ;
	}	
}

void Mileusnic06Spindle::extendInitStateFromProperties(SimTK::State& s) const
{	
	Super::extendInitStateFromProperties(s);
	setDynamicActivation(s,getDefaultActivation());
	setStaticActivation(s,getDefaultActivation());
	
	// also initiate the state variables that don't come from properties
	setTensionBag1(s,0.05);
	setTensionBag1Deriv(s,0.0);
	setTensionBag2(s,0.05);
	setTensionBag2Deriv(s,0.0);
	setTensionChain(s,0.05);
	setTensionChainDeriv(s,0.0);
}

void Mileusnic06Spindle::extendSetPropertiesFromState(const SimTK::State& s)
{	
	Super::extendSetPropertiesFromState(s);
	setDefaultActivation(getDynamicActivation(s)); // ignoring static activation
}

void Mileusnic06Spindle::extendConnectToModel(Model& aModel) 
{
	Super::extendConnectToModel(aModel);
	
	_model = &aModel;

	musclePtr = traversePathToComponent<Muscle>(getAbsolutePath().getParentPath());
}

//--------------------------------------------------------------------------
// GET & SET Properties
//--------------------------------------------------------------------------
void Mileusnic06Spindle::setDefaultActivation(double aDefaultActivation) {
	set_default_activation(aDefaultActivation);
}

void Mileusnic06Spindle::constructProperties()
{
	setAuthors("Sergio Verduzco");
	constructProperty_default_activation(0.05);
}

//--------------------------------------------------------------------------
// GET & SET States and their derivatives
//--------------------------------------------------------------------------

// dynamic_activation
double Mileusnic06Spindle::getDynamicActivation(const SimTK::State& s) const {
	return getStateVariableValue(s, STATE_dynamic_activation_NAME); 
}
void Mileusnic06Spindle::setDynamicActivation(SimTK::State& s, double Activation) const {
	setStateVariableValue(s, STATE_dynamic_activation_NAME, Activation); 
}
// static_activation
double Mileusnic06Spindle::getStaticActivation(const SimTK::State& s) const {
	return getStateVariableValue(s, STATE_static_activation_NAME); 
}
void Mileusnic06Spindle::setStaticActivation(SimTK::State& s, double Activation) const {
	setStateVariableValue(s, STATE_static_activation_NAME, Activation); 
}
// bag1 tension
double Mileusnic06Spindle::getTensionBag1(const SimTK::State& s) const {
	return getStateVariableValue(s, STATE_tension_bag1_NAME); 
}
void Mileusnic06Spindle::setTensionBag1(SimTK::State& s, double Tension) const {
	setStateVariableValue(s, STATE_tension_bag1_NAME, Tension); 
}
double Mileusnic06Spindle::getTensionBag1Deriv(const SimTK::State& s) const {
	return getStateVariableValue(s, STATE_tension_bag1_deriv_NAME); 
}
void Mileusnic06Spindle::setTensionBag1Deriv(SimTK::State& s, 
										double TensionDeriv) const {
	setStateVariableValue(s, STATE_tension_bag1_deriv_NAME, TensionDeriv);
}
// bag2 tension
double Mileusnic06Spindle::getTensionBag2(const SimTK::State& s) const {
	return getStateVariableValue(s, STATE_tension_bag2_NAME); 
}
void Mileusnic06Spindle::setTensionBag2(SimTK::State& s, double Tension) const {
	setStateVariableValue(s, STATE_tension_bag2_NAME, Tension); 
}
double Mileusnic06Spindle::getTensionBag2Deriv(const SimTK::State& s) const {
	return getStateVariableValue(s, STATE_tension_bag2_deriv_NAME); 
}
void Mileusnic06Spindle::setTensionBag2Deriv(SimTK::State& s, 
										double TensionDeriv) const {
	setStateVariableValue(s, STATE_tension_bag2_deriv_NAME, TensionDeriv);
}
// chain tension
double Mileusnic06Spindle::getTensionChain(const SimTK::State& s) const {
	return getStateVariableValue(s, STATE_tension_chain_NAME); 
}
void Mileusnic06Spindle::setTensionChain(SimTK::State& s, double Tension) const {
	setStateVariableValue(s, STATE_tension_chain_NAME, Tension); 
}
double Mileusnic06Spindle::getTensionChainDeriv(const SimTK::State& s) const {
	return getStateVariableValue(s, STATE_tension_chain_deriv_NAME); 
}
void Mileusnic06Spindle::setTensionChainDeriv(SimTK::State& s, 
										double TensionDeriv) const {
	setStateVariableValue(s, STATE_tension_chain_deriv_NAME, TensionDeriv);
}
// output variables
double Mileusnic06Spindle::getIaOutput(const SimTK::State& s) const {
	return getCacheVariableValue<double>(s, CACHE_primaryIa_NAME);
}
double& Mileusnic06Spindle::updIaOutput(const SimTK::State& s) const
{
    return updCacheVariableValue<double>(s, CACHE_primaryIa_NAME);
}
double Mileusnic06Spindle::getIIOutput(const SimTK::State& s) const {
	return getCacheVariableValue<double>(s, CACHE_secondaryII_NAME);
}
double& Mileusnic06Spindle::updIIOutput(const SimTK::State& s) const
{
    return updCacheVariableValue<double>(s, CACHE_secondaryII_NAME);
}

//=============================================================================
// COMPUTATION OF DERIVATIVES
//=============================================================================
 
void Mileusnic06Spindle::computeStateVariableDerivatives(const SimTK::State& s) const
{
	// vector of the derivatives to be returned
	SimTK::Vector derivs(getNumStateVariables(), SimTK::NaN);
	int nd = derivs.size();
	SimTK_ASSERT1(nd == 8, "Mileusnic06Spindle: Expected 8 state variables"
        " but encountered  %f.", nd);
	// The state variables corresponding to the entries in derivs are: 
	// derivs[0] --> Dynamic bag fiber activation
	// derivs[1] --> Static bag fiber activation
	// derivs[2] --> Dynamic bag fiber tension
	// derivs[3] --> Static bag fiber tension
	// derivs[4] --> Chain fiber tension
	// derivs[5] --> Dynamic bag fiber tension's derivative
	// derivs[6] --> Static bag fiber tension's derivative
	// derivs[7] --> Chain fiber tension's derivative
	
	// The fusimotor inputs come from an AfferentController object
	// that provides 3 controls to the muscle. 
	SimTK::VectorView ctrlVec = musclePtr->getControls(s);
	// ctrlVec[0] --> muscle excitation
	// ctrlVec[1] --> dynamic   fusimotor inputs
	// ctrlVec[2] --> static fusimotor inputs

	// Sometimes, controls can come in as NaNs - very concerning.
	// Let's assume it's supposed to be 0 for now.
	double gssq = 0.0;
	if (!isnan(ctrlVec[2]))
		gssq = (ctrlVec[2]) * (ctrlVec[2]);

	double gdsq = 0.0;
	if (!isnan(ctrlVec[1]))
		gdsq = (ctrlVec[1]) * (ctrlVec[1]);

	// Calculating activation derivatives
	double fdsq = (bag1.freq)*(bag1.freq); // I could store these two 
	double fssq = (bag2.freq)*(bag2.freq); // to speed things up
	double fcsq = (chain.freq)*(chain.freq);

	// Dynamic Activation (Equation 1)
	derivs[0] = ( (gdsq/(gdsq+fdsq)) - getDynamicActivation(s) )/bag1.tau;
	// Static Activation (Equation 1)
	derivs[1] = ( (gssq/(gssq+fssq)) - getStaticActivation(s) )/bag2.tau;
	// calculating chain fiber activation, which is instantaneous (Equation 1)
	double ch_act;
	ch_act = gssq / ( gssq + fcsq );
	
	// Calculating damping terms (beta)
	double beta_bag1, beta_bag2, beta_chain;
	// dynamic damping term. Equation 4 with beta_2 = 0.
	beta_bag1 = bag1.beta_0 + bag1.beta_1*getDynamicActivation(s);
	// static bag damping term. Equation 4 with beta_1 = 0.
	beta_bag2 = bag2.beta_0 + bag2.beta_2*getStaticActivation(s);
	// chain damping term. Equation 4 with beta_1 = 0.
	beta_chain = chain.beta_0 + chain.beta_2*ch_act;
	
	// calculating the force generator terms (Gamma)
	double Gamma_bag1, Gamma_bag2, Gamma_chain;
	// bag1 force generator term. Equation 5 with Gamma_2 = 0
	Gamma_bag1 = bag1.Gamma_1*getDynamicActivation(s);
	// bag2 force generator term. Equation 5 with Gamma_1 = 0
	Gamma_bag2 = bag2.Gamma_2*getStaticActivation(s);
	// chain force generator term. Equation 5 with Gamma_1 = 0
	Gamma_chain = chain.Gamma_2*ch_act;

	// The first derivative of the tension state variables
	// comes from another set of state variables
	derivs[2] = getTensionBag1Deriv(s);
	derivs[3] = getTensionBag2Deriv(s);
	derivs[4] = getTensionChainDeriv(s);
	
	// Now let's compute the second derivative of the tension
	// (Eq. 6), and also get the afferent potentials (Eqs. 7, 8)
	
	// equation 6 needs some muscle information
	double L0 = musclePtr->getOptimalFiberLength();
	double L = musclePtr->getFiberLength(s)/L0;
	double Lp = (musclePtr->getFiberVelocity(s))/L0;
	double Lpp = ((Millard12EqMuscleWithAfferents*)musclePtr)->
	             getLPFacceleration(s)/L0;
	double  C, term1, term2, T, Tp;	// auxiliary variables 
	
	// Tension 2nd derivative for bag1 
	C = (Lp>0.0) ? bag1.C_L : bag1.C_S;
	T = getTensionBag1(s);
	Tp = derivs[2];
	 
	term1 = C * beta_bag1 * sgn(Lp - (Tp/bag1.K_SR))
		  * std::pow(std::abs(Lp - (Tp/bag1.K_SR)), bag1.a)
		  * (L - bag1.L_0SR - (T/bag1.K_SR) - bag1.R);
	term2 = bag1.K_PR*(L - bag1.L_0SR - (T/bag1.K_SR) - bag1.L_0PR);
	
	derivs[5] = (bag1.K_SR/bag1.M) * ( term1 + term2 
	          + bag1.M*Lpp + Gamma_bag1 - T );
	 
	// afferent potential for bag1 (equation 7)
	double APbag1; 
	//APbag1 = bag1.G * max( (T/bag1.K_SR) - (bag1.L_NSR - bag1.L_0SR) , 0.0);
	APbag1 = bag1.G * ( (T/bag1.K_SR) - (bag1.L_NSR - bag1.L_0SR) );			
			
	// Tension 2nd derivative for bag2 
	C = (Lp>0.0) ? bag2.C_L : bag2.C_S;
	T = getTensionBag2(s);
	Tp = derivs[3];
	 
	term1 = C * beta_bag2 * sgn(Lp - (Tp/bag2.K_SR))
		  * std::pow(std::abs(Lp - (Tp/bag2.K_SR)), bag2.a)
		  * (L - bag2.L_0SR - (T/bag2.K_SR) - bag2.R);
	term2 = bag2.K_PR*(L - bag2.L_0SR - (T/bag2.K_SR) - bag2.L_0PR);
	
	derivs[6] = (bag2.K_SR/bag2.M) * ( term1 + term2 
	          + bag2.M*Lpp + Gamma_bag2 - T );
			
	// afferent potential for bag2 (equation 8 except G product)
	double APbag2;
	/*
	APbag2 = bag2.X * (bag2.L_sec/bag2.L_0SR)
	                * max((T/bag2.K_SR) - (bag2.L_NSR - bag2.L_0SR), 0.0)
		   + (1.0-bag2.X) * (bag2.L_sec/bag2.L_0PR)
					* max(L - (T/bag2.K_SR) - bag2.L_0SR - bag2.L_NPR, 0.0); 
	*/

	APbag2 = bag2.X * (bag2.L_sec/bag2.L_0SR)
	                * ((T/bag2.K_SR) - (bag2.L_NSR - bag2.L_0SR))
		   + (1.0-bag2.X) * (bag2.L_sec/bag2.L_0PR)
					* (L - (T/bag2.K_SR) - bag2.L_0SR - bag2.L_NPR); 
					
	// Tension 2nd derivative for the chain fiber
	C = (Lp>0.0) ? chain.C_L : chain.C_S;
	T = getTensionChain(s);
	Tp = derivs[4];
	 
	term1 = C * beta_chain * sgn(Lp - (Tp/chain.K_SR))
		  * std::pow(std::abs(Lp - (Tp/chain.K_SR)), chain.a)
		  * (L - chain.L_0SR - (T/chain.K_SR) - chain.R);	
	term2 = chain.K_PR*(L - chain.L_0SR - (T/chain.K_SR) - chain.L_0PR);
	
	derivs[7] = (chain.K_SR/chain.M) * ( term1 + term2 
	          + chain.M*Lpp + Gamma_chain - T );
			  
	// afferent potential for chain (equation 8 except G product)
	double APchain;
	/*
	APchain = chain.X * (chain.L_sec/chain.L_0SR)
	                  * max((T/chain.K_SR) - (chain.L_NSR - chain.L_0SR), 0.0)
		    + (1.0-chain.X) * (chain.L_sec/chain.L_0PR)
					  * max(L - (T/chain.K_SR) - chain.L_0SR - chain.L_NPR, 0.0);
	*/
	APchain = chain.X * (chain.L_sec/chain.L_0SR)
	                  * ((T/chain.K_SR) - (chain.L_NSR - chain.L_0SR))
		    + (1.0-chain.X) * (chain.L_sec/chain.L_0PR)
					  * (L - (T/chain.K_SR) - chain.L_0SR - chain.L_NPR);
					  
	// calculating the afferent firing
	double primary, secondary, pri_stat;
	pri_stat = bag2.G_pri*APbag2 + chain.G_pri*APchain;
	primary = max(APbag1, pri_stat) + S * min(APbag1, pri_stat);
	secondary = bag2.G_sec*APbag2 + chain.G_sec*APchain;
	
	// cache the output so it can be accessed
	double& iaout = updIaOutput(s);
    iaout = primary;
    markCacheVariableValid(s,CACHE_primaryIa_NAME);

	double& iiout = updIIOutput(s);
    iiout = secondary;
    markCacheVariableValid(s,CACHE_secondaryII_NAME);

	setStateVariableDerivativeValue(s, STATE_dynamic_activation_NAME, derivs[0]);
	setStateVariableDerivativeValue(s, STATE_static_activation_NAME, derivs[1]);
	setStateVariableDerivativeValue(s, STATE_tension_bag1_NAME, derivs[2]);
	setStateVariableDerivativeValue(s, STATE_tension_bag2_NAME, derivs[3]);
	setStateVariableDerivativeValue(s, STATE_tension_chain_NAME, derivs[4]);
	setStateVariableDerivativeValue(s, STATE_tension_bag1_deriv_NAME, derivs[5]);
	setStateVariableDerivativeValue(s, STATE_tension_bag2_deriv_NAME, derivs[6]);
	setStateVariableDerivativeValue(s, STATE_tension_chain_deriv_NAME, derivs[7]);
}

//--------------------------------------------------------------------------
// Method to set the initial conditions
//--------------------------------------------------------------------------
void Mileusnic06Spindle::
computeInitialSpindleEquilibrium(SimTK::State& s) const
{	
	// Calculating damping terms (beta)
	double beta_bag1, beta_bag2, beta_chain;
	// dynamic damping term. Equation 4 with beta_2 = 0.
	beta_bag1 = bag1.beta_0 + bag1.beta_1*getDynamicActivation(s);
	// static bag damping term. Equation 4 with beta_1 = 0.
	beta_bag2 = bag2.beta_0 + bag2.beta_2*getStaticActivation(s);
	// chain damping term. Equation 4 with beta_1 = 0.
	beta_chain = chain.beta_0 + chain.beta_2*getStaticActivation(s);
	
	// calculating the force generator terms (Gamma)
	double Gamma_bag1, Gamma_bag2, Gamma_chain;
	// bag1 force generator term. Equation 5 with Gamma_2 = 0
	Gamma_bag1 = bag1.Gamma_1*getDynamicActivation(s);
	// bag2 force generator term. Equation 5 with Gamma_1 = 0
	Gamma_bag2 = bag2.Gamma_2*getStaticActivation(s);
	// chain force generator term. Equation 5 with Gamma_1 = 0
	Gamma_chain = chain.Gamma_2*getStaticActivation(s);
	
	// normalized fiber length and velocity
	double L0 = musclePtr->getOptimalFiberLength();
	double L = musclePtr->getNormalizedFiberLength(s);
	double Lp = (musclePtr->getFiberVelocity(s))/L0;
	// clipping away large velocities 
	Lp = (Lp>15.0)? 15.0 : (Lp<-15.0)? -15.0 : Lp;
	
	//*******************************************************
	// Iteratively calculate LPR, dLPR, and tension
	double LPRb1, LPRb2, LPRc;
	double Tb1, Tb2, Tc, dTb1, dTb2, dTc;
	// Initial estimate of dLPR
	double dLPRb1, dLPRb2, dLPRc;
	dLPRb1 = bag1.K_SR * Lp / (bag1.K_SR + bag1.K_PR);
	dLPRb2 = bag2.K_SR * Lp / (bag2.K_SR + bag2.K_PR);
	dLPRc = chain.K_SR * Lp / (chain.K_SR + chain.K_PR);
	
	// calculate some terms that remain constant through
	// the iterations
	double Cb1, Cb2, Cc;
	Cb1 = (Lp>0.0) ? bag1.C_L : bag1.C_S;
	Cb2 = (Lp>0.0) ? bag2.C_L : bag2.C_S;
	Cc	= (Lp>0.0) ? chain.C_L : chain.C_S;
	
	double num_b1, num_b2, num_c;
	double den_b1, den_b2, den_c;
	num_b1 = bag1.K_SR*(L-bag1.L_0SR) + bag1.K_PR*bag1.L_0PR + Gamma_bag1;
	num_b2 = bag2.K_SR*(L-bag2.L_0SR) + bag2.K_PR*bag2.L_0PR + Gamma_bag2;
	num_c = chain.K_SR*(L-chain.L_0SR) + chain.K_PR*chain.L_0PR + Gamma_chain;
	den_b1 = bag1.K_SR + bag1.K_PR;
	den_b2 = bag2.K_SR + bag2.K_PR;
	den_c = chain.K_SR + chain.K_PR;
	
	double sig = (double)sgn(Lp);
	double rab1 = 1.0/bag1.a;
	double rab2 = 1.0/bag2.a;
	double rac = 1.0/chain.a;
	
	double work_pow; // working variable
	
	// iterations
	for( int i=0; i <= 4; i++)
	{
		// bag 1
		work_pow = beta_bag1*Cb1*sig*std::pow(std::abs(dLPRb1),bag1.a);
		LPRb1 = (num_b1 + work_pow) / (den_b1 + work_pow);
		Tb1 = bag1.K_SR*(L - LPRb1 - bag1.L_0SR);
		dLPRb1 = (Tb1 - bag1.K_PR*(LPRb1 - bag1.L_0PR) + Gamma_bag1)
		     / (beta_bag1*Cb1*(LPRb1 - bag1.R));
		dLPRb1 = sig*std::pow(std::abs(dLPRb1),rab1);
		
		// bag 2
		work_pow = beta_bag2*Cb2*sig*std::pow(std::abs(dLPRb2),bag2.a);
		LPRb2 = (num_b2 + work_pow) / (den_b2 + work_pow);
		Tb2 = bag2.K_SR*(L - LPRb2 - bag2.L_0SR);
		dLPRb2 = (Tb2 - bag2.K_PR*(LPRb2 - bag2.L_0PR) + Gamma_bag2)
		     / (beta_bag2*Cb2*(LPRb2 - bag2.R));
		dLPRb2 = sig*std::pow(std::abs(dLPRb2),rab2);
	
		// chain	
		work_pow = beta_chain*Cc*sig*std::pow(std::abs(dLPRc),chain.a);
		LPRc = (num_c + work_pow) / (den_c + work_pow);
		Tc = chain.K_SR*(L - LPRc - chain.L_0SR);
		dLPRc = (Tc - chain.K_PR*(LPRc - chain.L_0PR) + Gamma_chain)
		     / (beta_chain*Cc*(LPRc - chain.R));
		dLPRc = sig*std::pow(std::abs(dLPRc),rac);	
	}
	// TODO: check convergence
	
	// tension derivatives 
	dTb1 = bag1.K_SR * (Lp - dLPRb1);
	dTb2 = bag2.K_SR * (Lp - dLPRb2);
	dTc = chain.K_SR * (Lp - dLPRc);
	
	// set values
	setTensionBag1(s, Tb1);
	setTensionBag1Deriv(s, dTb1);
	setTensionBag2(s, Tb2);
	setTensionBag2Deriv(s, dTb2);
	setTensionChain(s, Tc);
	setTensionChainDeriv(s, dTc);
}
