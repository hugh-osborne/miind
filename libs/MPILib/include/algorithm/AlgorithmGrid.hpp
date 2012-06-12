// Copyright (c) 2005 - 2010 Marc de Kamps
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
//      If you use this software in work leading to a scientific publication, you should cite
//      the 'currently valid reference', which can be found at http://miind.sourceforge.net
#ifndef MPILIB_ALGORITHMS_ALGORITHMGRID_HPP_
#define MPILIB_ALGORITHMS_ALGORITHMGRID_HPP_

#include <vector>
#include <valarray>
#include <MPILib/include/BasicTypes.hpp>

namespace MPILib {

//! AlgorithmGrid
class AlgorithmGrid {
public:

	//! Create the state for a AlgorithmGrid with a maximum number of elements
	AlgorithmGrid(Number);

	//! copy constructor
	AlgorithmGrid(const AlgorithmGrid&);

	//! Create an Algorithmrid with just a state (usually a single number)
	AlgorithmGrid(const std::vector<double>&);

	//! construct an AlgorithmGrid from two a state and an interpretation
	AlgorithmGrid(const std::vector<double>&, const std::vector<double>&);

	AlgorithmGrid& operator=(const AlgorithmGrid&);

	std::vector<double> ToStateVector();
	std::vector<double> ToInterpretationVector();

private:

	template<class WeightValue> friend class AbstractAlgorithm;

	template<class Value>
	std::valarray<Value> ToValarray(const std::vector<double>& vector);

	template<class Value>
	std::vector<Value> ToVector(const std::valarray<Value>& array,
			Number number_to_be_copied);

	std::valarray<double>& ArrayState();
	std::valarray<double>& ArrayInterpretation();
	Number& StateSize();
	Number StateSize() const;

	void Resize(Number);

	//! allow iteration over internal values of the state
	const double* begin_state() const;
	const double* end_state() const;

	const double* begin_interpretation() const;
	const double* end_interpretation() const;

	Number _number_state; // the array_state is sometimes larger than the actual state
	std::valarray<double> _array_state;
	std::valarray<double> _array_interpretation;
};

} // end of namespace

#endif //MPILIB_ALGORITHMS_ALGORITHMGRID_HPP_ include guard
