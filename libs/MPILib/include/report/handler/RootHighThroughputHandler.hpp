// Copyright (c) 2005 - 2011 Marc de Kamps
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

#ifndef MPILIB_REPORT_HANDLER_ROOTHIGHTHROUGHPUTHANDLER_HPP_
#define MPILIB_REPORT_HANDLER_ROOTHIGHTHROUGHPUTHANDLER_HPP_

#include <string>
#include <vector>
#include <MPILib/include/report/handler/AbstractReportHandler.hpp>
#include <MPILib/include/BasicTypes.hpp>

class TTree;
class TFile;
class TBranch;
template <class Type> class TVectorT;

namespace MPILib {
namespace report {
namespace handler {
	//! This is a handler which organises data in time slices. Its memory use is constant over a simulation, unlike that
	//! of RootReportHandler, but is doesn't offer online visualisation.
	//!
	class RootHighThroughputHandler : public AbstractReportHandler
	{
	public:

		//! Standard constructor for client code
		RootHighThroughputHandler
		(
			const std::string&,	//! root file name
			bool,			//! write out state file 
			bool	reinstate_rate_graph = false //! backward compatibility option for an older ROOT layout
		);

		RootHighThroughputHandler(const RootHighThroughputHandler&);

		virtual ~RootHighThroughputHandler();


		//! Collects the Report of a DynamicNode for storage in the simulation file.
		virtual void writeReport(const Report&);


		virtual RootHighThroughputHandler* clone() const;

		virtual void initializeHandler
		(
				const NodeId&
		);


		virtual void detachHandler
		(
			const NodeId&
		);

	private:

		bool reinstateNodeGraphs(const char*);
	    void collectGraphInformation(std::vector<double>*,Number*,Number*);
		void storeRateGraphs(const std::vector<double>&, Number, Number);

		static Time					_t_start;
		static TTree*				_p_tree;
		static TFile*				_p_file;
		static TVectorT<double>*	_p_array;
		static bool					_reinstate_node_graphs;
		static bool					_is_recording;
		static bool					_is_first_time_slice_processed;
		static int					_nr_slice;
		static std::vector<double>		_vec_data;

	};

}// end namespace of handler
}// end namespace of report
}// end namespace of MPILib
#endif // include guard MPILIB_REPORT_HANDLER_ROOTHIGHTHROUGHPUTHANDLER_HPP_