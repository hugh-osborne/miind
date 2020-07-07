// Copyright (c) 2005 - 2009 Marc de Kamps, Dave Harrison
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
#include "pch.h"
#include "Interpolation.h"
#include "NumtoolsLibException.h"
#include <gsl/gsl_errno.h>

#include <iostream>

using namespace NumtoolsLib;

Interpolator::Interpolator
(
	const InterpType type,
	const std::valarray<double>& x, 
	const std::valarray<double>& y
) :
	_x(const_cast<std::valarray<double>&>(x)), // const cast is justfied by the fact that 
	_y(const_cast<std::valarray<double>&>(y)), // gsl_interp takes const pointers, valarray doesn't have a const method to access its data address
	_size(_x.size())
    {
    // Determine type and allocate an interpolator
    switch (type)
        {
        case INTERP_LINEAR:
            _interp = gsl_interp_alloc(gsl_interp_linear, _size);
            break;

        case INTERP_CSPLINE:
            _interp = gsl_interp_alloc(gsl_interp_cspline, _size);
            break;

        case INTERP_AKIMA:
            _interp = gsl_interp_alloc(gsl_interp_akima, _size);
            break;
        }

    // Create an accelerator object
    _acc = gsl_interp_accel_alloc();

    // Allocate an interpolator and initialise
    int gsl_error_code = gsl_interp_init(_interp, &_x[0], &_y[0], _size);
    if (gsl_error_code)
        throw NumtoolsException(gsl_strerror(gsl_error_code));
    }

double Interpolator::InterpValue(const double xi)
    {
    double yi;
    int gsl_error_code = gsl_interp_eval_e(_interp, &_x[0], &_y[0], xi, _acc,
            &yi);

    if (gsl_error_code)
        throw NumtoolsException(gsl_strerror(gsl_error_code));

    return yi;
    }

Interpolator::~Interpolator()
    {
    // Reset the accelerator
    if (_acc)
        gsl_interp_accel_reset(_acc);

    // Free up the resources for the interpolator
    if (_interp)
        gsl_interp_free(_interp);

    // Free up the resources for the accelerator
    if (_acc)
        gsl_interp_accel_free(_acc);
}

