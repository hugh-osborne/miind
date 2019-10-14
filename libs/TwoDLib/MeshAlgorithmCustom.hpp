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

#ifndef _CODE_LIBS_TWODLIB_MESHALGORITHMCUSTOM_INCLUDE_GUARD
#define _CODE_LIBS_TWODLIB_MESHALGORITHMCUSTOM_INCLUDE_GUARD

#include <string>
#include <vector>
#include <MPILib/include/AlgorithmInterface.hpp>
#include <MPILib/include/DelayedConnectionQueue.hpp>
#include <MPILib/include/CustomConnectionParameters.hpp>
#include "MasterOdeint.hpp"
#include "Ode2DSystemGroup.hpp"
#include "pugixml.hpp"
#include "display.hpp"
#include "DensityAlgorithmInterface.hpp"

namespace TwoDLib {

/**
 * \brief Mesh or 2D algorithm class.
 *
 * This class simulates the evolution of a neural population density function on a 2D grid.
 */

	template <class Solver=TwoDLib::MasterOdeint>
	class MeshAlgorithmCustom : public DensityAlgorithmInterface<CustomConnectionParameters>  {

	public:

		MeshAlgorithmCustom
		(
			const std::string&, 	                //!< model file name
			const std::vector<std::string>&,        //!< vector of transition matrix files
			MPILib::Time,                           //!< default time step for Master equation
			MPILib::Time tau_refractive = 0,        //!< absolute refractive period
			const string& ratemethod = ""           //!< firing rate computation; by default the mass flux across threshold
		);

		MeshAlgorithmCustom(const MeshAlgorithmCustom&);


		/**
		 * Cloning operation, to provide each DynamicNode with its own
		 * Algorithm instance. Clients use the naked pointer at their own risk.
		 */
		virtual MeshAlgorithmCustom* clone() const;

		/**
		 * Configure the Algorithm
		 * @param simParam The simulation parameter
		 */
		virtual void configure(const MPILib::SimulationRunParameter& simParam);


		/**
		 * The current time point
		 * @return The current time point
		 */
		virtual MPILib::Time getCurrentTime() const {return _t_cur;}


		/**
		 * The calculated rate of the node
		 * @return The rate of the node
		 */
	  virtual MPILib::Rate getCurrentRate() const {return _rate;}

		/**
		 * Stores the algorithm state in a Algorithm Grid
		 * @return The state of the algorithm
		 */
	  virtual MPILib::AlgorithmGrid getGrid(MPILib::NodeId, bool b_state = true) const;

		virtual void reportDensity(MPILib::Time t) const;

		/**
		 * Evolve the node state. In the default case it simply calls envolveNodeState
		 * without the NodeTypes. However if an algorithm needs the nodeTypes
		 * of the precursors overwrite this function.
		 * @param nodeVector Vector of the node States
		 * @param weightVector Vector of the weights of the nodes
		 * @param time Time point of the algorithm
		 * @param typeVector Vector of the NodeTypes of the precursors
		 */
		using MPILib::AlgorithmInterface<CustomConnectionParameters>::evolveNodeState;
		virtual void evolveNodeState(const std::vector<MPILib::Rate>& nodeVector,
				const std::vector<CustomConnectionParameters>& weightVector, MPILib::Time time);

		/**
		 * prepare the Evolve method
		 * @param nodeVector Vector of the node States
		 * @param weightVector Vector of the weights of the nodes
		 * @param typeVector Vector of the NodeTypes of the precursors
		 */
		virtual void prepareEvolve(const std::vector<MPILib::Rate>& nodeVector,
				const std::vector<CustomConnectionParameters>& weightVector,
				const std::vector<MPILib::NodeType>& typeVector);


		/**
		 * Provides a reference to the Mesh
		 */
		const Mesh& MeshReference() const { return _mesh_vec[0]; }

		/**
		 * By default, to find a matrix associated to an efficacy, MeshAlgorithm tests whether the efficacy quoted
		 * by a TransitionMatrix is within a certain tolerance of a given efficacy. If that tolerance is not appropriate,
		 * it can be adapted here
		 */
		void SetWeightTolerance(double tolerance){_tolerance = tolerance; }

		/**
		 * Obtain current weight tolerance
		 */
		double GetWeightTolerance() const { return _tolerance; }

		/**
		 * Initialize a given Mesh cell before simulation starts
		 */
		void InitializeDensity(MPILib::Index i, MPILib::Index j){_sys.Initialize(0,i,j);}

	  const Ode2DSystemGroup& Sys() const { return _sys; }


		//keep track of the node id
		virtual void assignNodeId( MPILib::NodeId );


		const std::vector<TwoDLib::Redistribution>& MapReversal() const { return  _vec_vec_rev[0]; }
		const std::vector<TwoDLib::Redistribution>& MapReset()    const { return  _vec_vec_res[0]; }

	private:


		// initialization routines
		pugi::xml_node                          CreateRootNode(const std::string&);
		std::vector<TwoDLib::Mesh>              CreateMeshObject();
		std::vector<TwoDLib::Redistribution>    Mapping(const std::string&);
		std::vector<TwoDLib::TransitionMatrix>  InitializeMatrices(const std::vector<std::string>&);
		void                                    FillMap(const std::vector<CustomConnectionParameters>& weightVector);

		double _tolerance;

		const std::string              _model_name;
		const std::vector<std::string> _mat_names;   // it is useful to store the names, but not the matrices, as they will be converted internally by the MasterOMP object
        const std::string              _rate_method;

		// report quantities
		MPILib::Time _h;
		MPILib::Rate _rate;
		MPILib::Time _t_cur;

		std::vector<MPILib::Time>    _vec_tau_refractive;

		// parsing auxilliaries
		pugi::xml_document _doc;
		pugi::xml_node     _root;

		MPILib::NodeId _node_id;

		// mesh and mappings
		std::vector<TwoDLib::Mesh> _mesh_vec;          // we have to create a vector, even if MeshAlgorithm manages only one Mesh, because
		                                               // there needs to be a place for Ode2DSystemGroup to store a vector of Mesh
		std::vector< std::vector<TwoDLib::Redistribution> > _vec_vec_rev;
		std::vector< std::vector<TwoDLib::Redistribution> > _vec_vec_res;

		// map incoming rates onto the order used by MasterOMP
		std::vector<MPILib::Index>              _vec_map;
		std::vector<std::vector<MPILib::DelayedConnectionQueue>> _vec_vec_delay_queues; // this is fed to the apply step of MasterOMP.

		MPILib::Time 						_dt;     // mesh time step
		MPILib::Time						_network_time_step; // network time step
		TwoDLib::Ode2DSystemGroup 			_sys;

		std::unique_ptr<Solver>	_p_master;
		MPILib::Number			_n_evolve;
		MPILib::Number			_n_steps;

		const vector<double>& (TwoDLib::Ode2DSystemGroup::*_sysfunction) () const;
	};
}

#endif // MeshAlgorithm
