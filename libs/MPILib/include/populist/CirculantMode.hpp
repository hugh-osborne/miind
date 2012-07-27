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
#ifndef MPILIB_POPULIST_CIRCULANTMODE_HPP_
#define MPILIB_POPULIST_CIRCULANTMODE_HPP_

namespace MPILib {
namespace populist {

//! In the Configure method of AbstractCirculantSolver an InputSetParameter reference is passed in. This contains
//! both an integer interpretation of the current potential jump in terms of number of bins (e.g. _H_exc) as well
//! as a floating point. DiffusionZeroLeakEquations will work with an integer version, whilst now for example
//! SingleInputZeroLeakEquations allow probability transport from one bin to a point between two bins, which requires
//! a floating point interpretation of the step size. This choice must be taken by a ZeroLeak developer. This developer must also
//! ensure that the CirculantMode and the NonCirculatMode are used consistently. Therefore also AbstractNonCirculantSolver uses this enum.
enum CirculantMode {
	FLOATING_POINT, INTEGER
};
} //end of namespace populist
} // end of namespace MPILib

#endif // include guard MPILIB_POPULIST_CIRCULANTMODE_HPP_
