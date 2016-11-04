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

#ifndef _CODE_LIBS_TWODLIBLIB_MESHALGORITHMCODE_INCLUDE_GUARD
#define _CODE_LIBS_TWODLIBLIB_MESHALGORITHMCODE_INCLUDE_GUARD

#include <boost/filesystem.hpp>
#include <MPILib/include/utilities/Log.hpp>
#include <fstream>
#include <iostream>
#include "MeshAlgorithm.hpp"
#include "Stat.hpp"
#include "TwoDLibException.hpp"

namespace {

	// predicate that helps to locate "Mapping" nodes in an xml structure
	class Pred {
	public:
		Pred(const std::string& type):_type(type){}

		bool operator() (pugi::xml_node node){

			return (std::string(node.name()) == "Mapping" && std::string(node.attribute("type").value()) == _type ) ? true : false;
		}

	private:
		std::string _type;
	};
}

namespace TwoDLib {

	template <class WeightValue>
	pugi::xml_node MeshAlgorithm<WeightValue>::CreateRootNode(const string& model_name){

		// document
		pugi::xml_parse_result result = _doc.load_file(model_name.c_str());
		pugi::xml_node  root = _doc.first_child();

		if (result.status != pugi::status_ok)
		  throw TwoDLib::TwoDLibException("Can't open .model file.");
		return root;
	}

	template <class WeightValue>
	Mesh MeshAlgorithm<WeightValue>::CreateMeshObject(){
		// mesh
		pugi::xml_node mesh_node = _root.first_child();

		if (mesh_node.name() != std::string("Mesh") )
		  throw TwoDLib::TwoDLibException("Couldn't find mesh node in model file");
		std::ostringstream ostmesh;
		mesh_node.print(ostmesh);
		std::istringstream istmesh(ostmesh.str());

		TwoDLib::Mesh mesh(istmesh);

		// MatrixGenerator should already have inserted the stationary bin and there is no need
		// to reexamine the stat file

		return mesh;
	}


	template <class WeightValue>
	std::vector<TwoDLib::Redistribution> MeshAlgorithm<WeightValue>::Mapping(const string& type)
	{
		Pred pred(type);
		pugi::xml_node rev_node = _root.find_child(pred);

		if (rev_node.name() != std::string("Mapping") ||
		    rev_node.attribute("type").value() != type)
			throw TwoDLibException("Couldn't find mapping in model file");

		std::ostringstream ostrev;
		rev_node.print(ostrev);
		std::istringstream istrev(ostrev.str());
		vector<TwoDLib::Redistribution> vec_rev = TwoDLib::ReMapping(istrev);
		return vec_rev;
	}

	template <class WeightValue>
	MeshAlgorithm<WeightValue>::MeshAlgorithm
	(
		const std::string& model_name,
		const std::vector<std::string>& mat_names,
		MPILib::Time h):
	_tolerance(1e-7),
	_model_name(model_name),
	_mat_names(mat_names),
	_h(h),
	_rate(0.0),
	_t_cur(0.0),
	_root(CreateRootNode(model_name)),
	_mesh(CreateMeshObject()),
	_vec_rev(this->Mapping("Reversal")),
	_vec_res(this->Mapping("Reset")),
	_vec_map(0),
	_dt(_mesh.TimeStep()),
	_sys(_mesh,_vec_rev,_vec_res),
	_n_evolve(0),
	_n_steps(0)
	// master parameter can only be calculated on configuration
	{
		// default initialization is (0,0); if there is no strip 0, it's down to the user
		if (_mesh.NrCellsInStrip(0) > 0 )
			_sys.Initialize(0,0);
	}

	template <class WeightValue>
	MeshAlgorithm<WeightValue>::MeshAlgorithm(const MeshAlgorithm<WeightValue>& rhs):
	_tolerance(rhs._tolerance),
	_model_name(rhs._model_name),
	_mat_names(rhs._mat_names),
	_h(rhs._h),
	_rate(rhs._rate),
	_t_cur(rhs._t_cur),
	_mesh(rhs._mesh),
	_vec_rev(rhs._vec_rev),
	_vec_res(rhs._vec_res),
	_vec_map(0),
	_dt(_mesh.TimeStep()),
	_sys(_mesh,_vec_rev,_vec_res),
	_n_evolve(0),
	_n_steps(0)
	// master parameter can only be calculated on configuration
	{
		// default initialization is (0,0); if there is no strip 0, it's down to the user
		if (_mesh.NrCellsInStrip(0) > 0 )
			_sys.Initialize(0,0);
	}

	template <class WeightValue>
	std::vector<TwoDLib::TransitionMatrix> MeshAlgorithm<WeightValue>::InitializeMatrices(const std::vector<std::string>& mat_names)
	{
		std::vector<TwoDLib::TransitionMatrix> vec_mat;

		for (const auto& name: mat_names)
			vec_mat.push_back(TransitionMatrix(name));
		return vec_mat;
	}

	template <class WeightValue>
	MeshAlgorithm<WeightValue>* MeshAlgorithm<WeightValue>::clone() const
	{
	  return new MeshAlgorithm<WeightValue>(*this);
	}

	template <class WeightValue>
	void MeshAlgorithm<WeightValue>::configure(const MPILib::SimulationRunParameter& par_run)
	{
		_t_cur = par_run.getTBegin();
		MPILib::Time t_step     = par_run.getTStep();

		// the integration time step, stored in the MasterParameter, is gauged with respect to the
		// network time step.
		MPILib::Number n_ode = static_cast<MPILib::Number>(std::floor(t_step/_h));
		MasterParameter par((n_ode > 1) ? n_ode : 1 );

		// vec_mat will go out of scope; MasterOMP will convert the matrices
		// internally and we don't want to keep two versions.
		std::vector<TransitionMatrix> vec_mat = InitializeMatrices(_mat_names);
 		std::unique_ptr<TwoDLib::MasterOMP> p_master(new MasterOMP(_sys,vec_mat, par));

		_p_master = std::move(p_master);

		// at this stage initialization must have taken place, either by default in (0,0),
		// or by the user calling Initialize if there is no strip 0

		double sum = _sys.P();
		if (sum == 0.)
			throw TwoDLib::TwoDLibException("No initialization of the mass array has taken place. Call Initialize before configure.");
	}

	template <class WeightValue>
	MPILib::AlgorithmGrid MeshAlgorithm<WeightValue>::getGrid(MPILib::NodeId id) const
	{
		// An empty grid will lead to crashes
		vector<double> array_interpretation {0.};
		vector<double> array_state {0.};

		std::ostringstream ost;
		ost << id << "_" << _t_cur;
		string fn("mesh_" + ost.str());

		if (!boost::filesystem::exists(_model_name + "_mesh" ) )
			boost::filesystem::create_directory(_model_name + "_mesh");

		std::ofstream ofst(_model_name + "_mesh/" + fn);
		_sys.Dump(ofst);

		return MPILib::AlgorithmGrid(array_state,array_interpretation);
	}

	template <class WeightValue>
	void MeshAlgorithm<WeightValue>::evolveNodeState
	(
		const std::vector<MPILib::Rate>& nodeVector,
		const std::vector<WeightValue>& weightVector,
		MPILib::Time time,
		const std::vector<MPILib::NodeType>& typeVector
	)
	{

	  // The network time step must be an integer multiple of the network time step; in principle
	  // we would expect this multiple to be one, but perhaps there are reasons to allow a population
	  // have a finer time resolution than others, so we allow larger multiples but write a warning in the log file.
		// determine number of steps and fix at the first step.
		if (_n_steps == 0){
		  // since n_steps == 0, time is the network time step
			double n = (time - _t_cur)/_dt;
		    
			_n_steps = static_cast<MPILib::Number>(round(n));
			if (_n_steps == 0){

			  throw TwoDLibException("Network time step is smaller than this grid's time step.");
			}
			if (fabs(_n_steps - n) > 1e-6){
			  throw TwoDLibException("Mismatch of mesh time step and network time step. Network time step should be a multiple (mostly one) of network time step");
			}
			if (_n_steps > 1)
			  LOG(MPILib::utilities::logWARNING)<< "Mesh runs at a time step which is a multiple of the network time step. Is this intended?";
			else
			  ; // else is fine
		}

		double evol;
	    // mass rotation
	    for (MPILib::Index i = 0; i < _n_steps; i++){
	      _sys.Evolve();
          _sys.RemapReversal();
          double reversal = _sys.P();
	    }

	    // master equation
	    _p_master->Apply(_n_steps*_dt,_vec_mapped_rates);
        _sys.RedistributeProbability();


 	    _t_cur += _n_steps*_dt;
 	    _rate = _sys.F();
 	    _n_evolve++;
	}

	template <class WeightValue>
	void MeshAlgorithm<WeightValue>::FillMap(const std::vector<WeightValue>& vec_weights)
	{
		// this function will only be called once;

 		for(const auto& weight: vec_weights){
			for (unsigned int i = 0; i < _p_master->NrMatrix(); i++){
				if ( fabs( _p_master->Efficacy(i) - weight._efficacy) < _tolerance ){
					if (std::find(_vec_map.begin(),_vec_map.end(),i) == _vec_map.end())
						_vec_map.push_back(i);
					else{
						throw TwoDLib::TwoDLibException("Weight in FillMap not unique");
					}
				}
			}
		}

		if (_vec_map.size() != vec_weights.size()){
			throw TwoDLib::TwoDLibException("FillMap map size differs from number of weights");
		}
		_vec_mapped_rates = std::vector<MPILib::Rate>(_vec_map.size(),0);
	}

	template <class WeightValue>
	void MeshAlgorithm<WeightValue>::prepareEvolve
	(
		const std::vector<MPILib::Rate>& nodeVector,
		const std::vector<WeightValue>& weightVector,
		const std::vector<MPILib::NodeType>& typeVector
	)
	{
		if (weightVector.size() != _p_master->NrMatrix())
			throw TwoDLib::TwoDLibException("There must be as many connections to a population as there are matrices.");
		// We assume that: 1) the ordering of weights during an evolve never changes; 2) transition matrices can be uniquely defined by their efficacies

		if (_vec_map.size() == 0)
			FillMap(weightVector);

		for (unsigned int i = 0; i < nodeVector.size(); i++){
			_vec_mapped_rates[_vec_map[i]] = nodeVector[i];
		}
	}
}

#endif // include guard
