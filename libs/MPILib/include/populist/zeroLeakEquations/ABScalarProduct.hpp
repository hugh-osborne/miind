// Copyright (c) 2005 - 2012 Marc de Kamps
//						2012 David-Matthias Sichau
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
#ifndef MPILIB_POPULIST_ZEROLEAKEQUATIONS_ABSCALARPRODUCT_HPP_
#define MPILIB_POPULIST_ZEROLEAKEQUATIONS_ABSCALARPRODUCT_HPP_

#include <MPILib/config.hpp>
#include <MPILib/include/populist/zeroLeakEquations/ABQStruct.hpp>
#include <MPILib/include/populist/OrnsteinUhlenbeckConnection.hpp>
#include <MPILib/include/utilities/Exception.hpp>
#include <vector>

namespace MPILib {
namespace populist {
namespace zeroLeakEquations {
/**
 * @brief calculates the AB Scalar Product
 */
class ABScalarProduct {
public:
	/**
	 * Evaluate the inner product over connections
	 * @param nodeVector Vector of the Rates of the precursor
	 * @param weightVector Vector of the weights of the connections
	 * @param time The current time
	 * @return A ABQ struct
	 */
	ABQStruct Evaluate(const std::vector<Rate>& nodeVector,
			const std::vector<OrnsteinUhlenbeckConnection>& weightVector,
			Time time) const {
		// for now a comes from the first population and b from the second. This will need to change.
		ABQStruct par_ret;

		if ((*nodeVector.begin()) == 6.0) {
			par_ret._a = 6.91423056;
			par_ret._b = 0.13299526;
		} else if ((*nodeVector.begin()) == 8.0) {
			par_ret._a = 129.43365395;
			par_ret._b = 0.08430153;
		} else
			throw utilities::Exception(
					"Input rate cannot be decoded by ABScalarProduct");

		return par_ret;
	}
};
} /* namespace zeroLeakEquations */
} /* namespace populist */
} /* namespace MPILib */

#endif // include guard MPILIB_POPULIST_ZEROLEAKEQUATIONS_ABSCALARPRODUCT_HPP_