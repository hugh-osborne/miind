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
//      If you use this software in work leading to a scientific publication, you should include a reference there to
//      the 'currently valid reference', which can be found at http://miind.sourceforge.net

#ifndef _CODE_2DLIB_ODE2DSYSTEMGROUP_INCLUDE_GUARD
#define _CODE_2DLIB_ODE2DSYSTEMGROUP_INCLUDE_GUARD

#include <iostream>
#include <numeric>
#include <vector>
#include "MPILib/include/TypeDefinitions.hpp"
#include "Mesh.hpp"
#include "modulo.hpp"
#include "Redistribution.hpp"
#include "CSRMatrix.hpp"
#include "MPILib/include/ProbabilityQueue.hpp"
#include "MPILib/include/RefractoryQueue.hpp"

namespace TwoDLib {

	class TransitionMatrix;


	  /**
	   * \brief Responsible for representing the density and deterministic evolution.
	   *
	   * This create an object that is similar in functionality to an Ode2DSystem, but for a group
	   * of meshes, so that the densities are stored in a single array, with the aim of vectorizing an
	   * Algorithm. You are strongly encouraged to test the Meshes and Mappings in a stand alone MeshAlgorithm simulation
	   * as this algorithm has more extensive checking for consistency between mesh and mappings.
	   */

	class Ode2DSystemGroup {
	public:

		//! Standard Constructor
		Ode2DSystemGroup
		(
			const std::vector<Mesh>&, 					       //!< A series of Mesh in the Python convention. Most models require a reversal bin that is not part of the grid. In that case it must be inserted into the Mesh by calling Mesh::InsertStationary. It is legal not to define an extra reversal bin, and use one of the existing Mesh cells at such, but in that case Cell (0,0) will not exist.
			const std::vector< std::vector<Redistribution> >&, //!< A series of mappings from strip end to reversal bin
			const std::vector< std::vector<Redistribution> >&,  //!< A series of mappings from threshold to reset bin
			const std::vector<MPILib::Time>&
		);

		Ode2DSystemGroup
		(
			const std::vector<Mesh>&, 					       //!< A series of Mesh in the Python convention. Most models require a reversal bin that is not part of the grid. In that case it must be inserted into the Mesh by calling Mesh::InsertStationary. It is legal not to define an extra reversal bin, and use one of the existing Mesh cells at such, but in that case Cell (0,0) will not exist.
			const std::vector< std::vector<Redistribution> >&, //!< A series of mappings from strip end to reversal bin
			const std::vector< std::vector<Redistribution> >&  //!< A series of mappings from threshold to reset bin
		);


		//! Place all initial density in a cell with coordinates (m,i,j)
		void Initialize(MPILib::Index, MPILib::Index,MPILib::Index);

		//! Map coordinates to a position in the density array. Map(0,0) may be defined or not, this depends on the Mesh and whether Mesh::InsertStationary
		//! was called. The safe way to use Map is in a loop that uses Mesh::NrCellsInStrip and Mesh::NrQuadrilateralStrips.
		MPILib::Index Map(
				unsigned int kmesh,   //! mesh number
				unsigned int istrip,  //! strip number
				unsigned int jcell    //! cell number
				) const {return _map[kmesh][istrip][jcell];}

		Coordinates toCoords(unsigned int index) const { return _index_to_coords[index]; }

		//! Takes an index for the mass array and converts to the mapped version. This is used in the Master equation solvers
		MPILib::Index Map(MPILib::Index i) const{ return _linear_map[i]; }

		//! Shift the density
		void Evolve();
		void Evolve(std::vector<MPILib::Index>& meshes);

		void EvolveWithoutMeshUpdate();

		//! Dump the current density profile (0), or the mass profile (1) to an output stream
		void Dump(const std::vector<std::ostream*>&, int mode = 0) const;

		//! Dump a single mesh file
		void DumpSingleMesh(std::ostream* vecost, unsigned int m, int mode = 0) const;

		//! Remap probability that has run from the end of a strip. Run this after evolution
		void RemapReversal();

		//! Redistribute probability that has moved through threshold. Run this after the Master equation
		void RedistributeProbability();

		void RedistributeProbability(MPILib::Number);

		void RedistributeIndividual(unsigned int threshold, unsigned int reset, unsigned int w_reset, MPILib::Number steps);

		//! Return the instantaneous firing rate
		const vector<MPILib::Rate>& F() const {return _fs;}

		//! total probability mass in the system, should not be too far away from 1.0
		MPILib::Mass P() const { return std::accumulate(_vec_mass.begin(),_vec_mass.end(),0.0); }

		//! average membrane potential, used for non-threshold crossing models such as Fitzhugh-Nagumo
		const vector<MPILib::Potential>& AvgV() const ;

		//! allow direct inspection of the mass array; client must still convert this to a density
		const vector<MPILib::Mass>& Mass() const { return _vec_mass; }

		//! allow inspection of the Mesh objects
		const std::vector<Mesh>& MeshObjects() const { return _mesh_list; }

		//! See what part of the mass array each Mesh is responsible for: the offsets are given as a function of mesh index
		const std::vector<MPILib::Index>& Offsets() const {return _vec_mesh_offset; }

		//! Provide access to the mass array
	    vector<MPILib::Mass>& Mass() { return _vec_mass; }

		const std::vector<double>& Vs() const { return _vec_vs; }

		void setIndividuals(vector<unsigned int> individuals) {_individuals = individuals; }

	    //! Provide read access to the Reversal map
		const std::vector<std::vector<Redistribution> >& MapReversal() const { return  _vec_reversal;}

	    //! Provide read access to the Reversal map
		const std::vector<std::vector<Redistribution> >& MapReset()    const { return  _vec_reset;}

		const std::vector<MPILib::Time> Tau_ref() const { return _vec_tau_refractive; }

		// Initialize the refractory reset queue with the network's time step
		void    InitializeResetRefractive(MPILib::Time network_time_step);

		friend class Master;
	    friend class MasterOMP;
	    friend class MasterOdeint;

	    friend void CheckSystem(const Ode2DSystemGroup&, const TransitionMatrix&, const std::vector<Redistribution>&, const std::vector<Redistribution>&, double);


	private:

		//! Copy constructor
		Ode2DSystemGroup(const Ode2DSystemGroup&);

		//! Implement the remapping of probability mass that hits the end of a strip
		class Reversal {
		public:

			Reversal(Ode2DSystemGroup& sys, vector<MPILib::Mass>& vec_mass, MPILib::Index m):_sys(sys),_vec_mass(vec_mass),_m(m){}

			void operator()(const Redistribution& map){
				_vec_mass[_sys.Map(_m,map._to[0],map._to[1])] += _vec_mass[_sys.Map(_m,map._from[0],map._from[1])];
				_vec_mass[_sys.Map(_m,map._from[0],map._from[1])] = 0;
			}

		private:
			Ode2DSystemGroup& 		_sys;
			vector<MPILib::Mass>&   _vec_mass;
			MPILib::Index     		_m;
		};

		//! Implement the remapping of probability mass that hits threshold
		class Reset {
		public:
			Reset(Ode2DSystemGroup& sys, vector<MPILib::Mass>& vec_mass,MPILib::Index m):_sys(sys),_vec_mass(vec_mass),_m(m){}

			void operator()(const Redistribution& map){
				MPILib::Mass from =  map._alpha*_vec_mass[_sys.Map(_m,map._from[0],map._from[1])];
				_vec_mass[_sys.Map(_m,map._to[0],map._to[1])] += from;
				_sys._fs[_m] += from;
			}

		private:

			Ode2DSystemGroup&		_sys;
			vector<MPILib::Mass>&  	_vec_mass;
			MPILib::Index       	_m;
		};

		class ResetRefractive {
		public:

			ResetRefractive
			(
				Ode2DSystemGroup&                  sys,
				vector<double>&               vec_mass,
				MPILib::Time                  network_time_step,
				MPILib::Time                  tau_refractive,
				const vector<Redistribution>& vec_reset,
				MPILib::Index m
			):
			_t_step(network_time_step),
			_tau_refractive(tau_refractive),
			_vec_reset(vec_reset),
			_vec_queue(vec_reset.size(),MPILib::RefractoryQueue(network_time_step,tau_refractive)),
			_sys(sys),
			_vec_mass(vec_mass),
			_m(m)
			{
			}

			void operator()(const Redistribution& map)
			{
				double from =  map._alpha*_vec_mass[_sys.Map(_m,map._from[0],map._from[1])];
				_sys._fs[_m] += from;

				MPILib::Index i = &map - &_vec_reset[0];
				from = _vec_queue[i].updateQueue(from);

				_vec_mass[_sys.Map(_m,map._to[0],map._to[1])] += from;
			}

			double getTotalProbInRefract() const{
				double total = 0.0;
				for (MPILib::RefractoryQueue q : _vec_queue){
					total += q.getTotalMass();
				}
				return total;
			}

		private:

			MPILib::Index       	_m;

			MPILib::Time                               _t_step;
			MPILib::Time                               _tau_refractive;

			const vector<Redistribution>&              _vec_reset;
			vector<MPILib::RefractoryQueue>  					_vec_queue;

			Ode2DSystemGroup&	                           _sys;
			vector<double>&                            _vec_mass;


		};

		//! Implement cleaning of the probability that was at threshold. TODO: this is mildly inefficient,
		//! in a serial implementation, but seems simpler in threads

		class Clean {
		public:

			Clean(Ode2DSystemGroup& sys, vector<MPILib::Mass>& vec_mass, MPILib::Index m):_sys(sys),_vec_mass(vec_mass),_m(m){
			}

			void operator()(const Redistribution& map){
				_vec_mass[_sys.Map(_m, map._from[0],map._from[1])] = 0.;
			}

		private:

			Ode2DSystemGroup&		_sys;
			vector<MPILib::Mass>&	_vec_mass;
			MPILib::Index       	_m;
		};

		vector<MPILib::Index> InitializeLength(const Mesh&) const;
		std::vector<std::vector<MPILib::Index> > InitializeLengths(const std::vector<Mesh>&);
		std::vector<MPILib::Index> InitializeCumulative(const Mesh&) const;
		std::vector<std::vector<MPILib::Index> > InitializeCumulatives(const std::vector<Mesh>&);
		std::vector<MPILib::Number> MeshOffset(const std::vector<Mesh>&) const;
		std::vector<double> MeshVs(const std::vector<Mesh>&) const;

		bool				  CheckConsistency() const;
		std::vector<Reset>    InitializeReset();
		std::vector<Reversal> InitializeReversal();
		std::vector<ResetRefractive> InitializeResetRefractiveInternal(MPILib::Time network_time_step);
		std::vector<Clean>    InitializeClean();

		vector<MPILib::Potential>  	InitializeArea(const std::vector<Mesh>&) const;
		vector<MPILib::Mass>        InitializeMass() const;

		std::vector< std::vector< std::vector<MPILib::Index> > > InitializeMap() const;
		std::vector<Coordinates> InitializeCoords() const;
		std::vector< MPILib::Index> InitializeLinearMap();
		std::vector< MPILib::Index> InitializeWorkingIndex();
		void 									UpdateMap(std::vector<MPILib::Index>& meshes);
		void                  UpdateMap();

		const std::vector<Mesh>&    _mesh_list;

		std::vector<MPILib::Number>              _vec_mesh_offset;
		std::vector<std::vector<MPILib::Index> > _vec_length;
		std::vector<std::vector<MPILib::Index> > _vec_cumulative;
		std::vector<double>					_vec_vs;

		std::vector<MPILib::Time>    _vec_tau_refractive;
public:
		vector<MPILib::Mass>	     	_vec_mass;
		vector<unsigned int> _individuals;
private:
		vector<MPILib::Potential>		_vec_area;

		unsigned int					_t;
		std::vector<MPILib::Mass>		_fs;
		std::vector<MPILib::Potential>	_avs;

		std::vector< std::vector<std::vector<MPILib::Index> > > _map;
		std::vector< MPILib::Index> _linear_map;
	  std::vector< Coordinates > _index_to_coords;

		std::vector<std::vector<Redistribution> > _vec_reversal;
		std::vector<std::vector<Redistribution> > _vec_reset;

		std::vector<Reversal> _reversal;
		std::vector<ResetRefractive> _reset_refractive;
		std::vector<Reset>    _reset;
		std::vector<Clean>    _clean;
	};
}
#endif // include guard
