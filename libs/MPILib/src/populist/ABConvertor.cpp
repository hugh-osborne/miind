#include "ABConvertor.h"
#include "AbstractCirculantSolver.h"
#include "PopulistSpecificParameter.h"


using namespace PopulistLib;


ABConvertor::ABConvertor
(
	VALUE_REF_INIT
	SpecialBins&,
	PopulationParameter&		par_pop,
	PopulistSpecificParameter&	par_specific,
	Potential&					delta_v,	
	Number&						n_current_bins
):
VALUE_MEMBER_INIT 
_p_specific	(&par_specific),
_p_pop		(&par_pop),
_p_n_bins	(&n_current_bins),
_p_delta_v	(&delta_v)
{
}

const PopulistSpecificParameter&
	ABConvertor::PopSpecific() const
{
	return *_p_specific;
}

const OneDMInputSetParameter&
	ABConvertor::InputSet() const
{
	return _param_input;
}

void ABConvertor::SortConnectionvector
(
	predecessor_iterator	iter_begin,
	predecessor_iterator	iter_end
)
{
	_param_input._par_input =
			_scalar_product.Evaluate
			(
				iter_begin, 
				iter_end,
				_p_pop->_tau	
			);
		_param_input._par_input._q = _param_onedm._par_adapt._q;
}

void ABConvertor::AdaptParameters
(	

)
{

	RecalculateSolverParameters();
}

void ABConvertor::RecalculateSolverParameters()
{
	_param_input._n_current_bins = *_p_n_bins;
	_param_input._n_max_bins     = _p_specific->MaxNumGridPoints();

	// current expansion factor is current number of bins
	// divided by number of initial bins

	double f = static_cast<double>(_param_input._n_current_bins)/static_cast<double>(_p_specific->NrGridInitial());
	_param_input._q_expanded = f*_param_input._par_input._q;
	_param_input._t_since_rebinning = _p_pop->_tau*log(f);
	_param_input._g_max = _p_pop->_theta;
	_param_input._tau   = _p_pop->_tau;
}