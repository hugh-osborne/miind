// Copyright (c) 2005 - 2015 Marc de Kamps
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation
//      and/or other materials provided with the distribution.
//    * Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software
//      without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#ifndef _CODE_LIBS_TWODLIB_MASTERGRID_INCLUDE_GUARD
#define _CODE_LIBS_TWODLIB_MASTERGRID_INCLUDE_GUARD

#include <string>
#include <boost/numeric/odeint.hpp>
#include "CSRMatrix.hpp"
#include "TransitionMatrix.hpp"
#include "Ode2DSystemGroup.hpp"
#include "MasterParameter.hpp"
#include <unordered_set>

namespace TwoDLib {

	//! OpenMP version of a forward Euler integration of the Master equation

	class MasterGrid {
	public:

		MasterGrid
		(
			Ode2DSystemGroup&,
			double
		);

		void MVGrid(
			vector<double>&       dydt,
			const vector<double>& vec_mass,
			double                rate,
		  double stays,
		  double goes,
		  int offset_1,
			int offset_2) const;

		void MVGridWithEfficacy(
			vector<double>&       dydt,
			const vector<double>& vec_mass,
			double                rate,
			unsigned int          efficiacy_index) const;

		void CalculateStaticEfficiacies(vector<double>& efficacy_map);
		void CalculateStaticEfficiaciesForConductance(vector<double>& efficacy_map, vector<double>& rest_v);

		void CalculateWindows(double t_step, const vector<double>& rates, vector<double>& efficacy_map);
		void Convolve(double t_step, const vector<double>& rates, vector<double>& efficacy_map, std::unordered_set<unsigned int>& cell_indices);

		void Apply(double t_step, const vector<double>& rates, vector<double>& efficacy_map);

		void operator()(const vector<double>&, vector<double>&, const double t = 0);

	private:

		MasterGrid& operator=(const MasterGrid&);

		Ode2DSystemGroup& _sys;

		double _cell_width;

		vector<double>			_dydt;
		vector<vector<double>>		_stays;
		vector<vector<double>>		_goes;
		vector<vector<int>>				_offset1;
		vector<vector<int>>				_offset2;

		unsigned int _window_size;
		vector<double> _dydt_window;
		vector<double> _mass_window;
		vector<unsigned int> add_inds;
    vector<unsigned int> subtract_inds;

		const vector<double>* _p_vec_eff;
		const vector<double>* _p_vec_rates;
	};
}

#endif // include guard
