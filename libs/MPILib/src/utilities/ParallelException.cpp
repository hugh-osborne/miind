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
//#include <MPILib/config.hpp>
#include <MPILib/include/utilities/ParallelException.hpp>
#include <MPILib/include/utilities/MPIProxy.hpp>
#include <sstream>

using namespace MPILib::utilities;

ParallelException::ParallelException(const char* message) :
		Exception(message) {
	std::stringstream sstream;
	sstream << std::endl << "Parallel Exception on processor: "
			<< MPIProxy().getRank() << " from: "
			<< MPIProxy().getSize()
			<< " with error message: " << msg_ << std::endl;
	msg_ = sstream.str();
}

ParallelException::ParallelException(const std::string& message) :
		Exception(message) {
	std::stringstream sstream;
	sstream << std::endl << "Parallel Exception on processor: "
			<< MPIProxy().getRank() << " from: "
			<< MPIProxy().getSize()
			<< " with error message: " << msg_ << std::endl;
	msg_ = sstream.str();
}

ParallelException::~ParallelException() throw(){
}

const char* ParallelException::what() const throw () {
	return msg_.c_str();
}
