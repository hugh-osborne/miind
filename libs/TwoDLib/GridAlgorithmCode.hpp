#ifndef _CODE_LIBS_TWODLIBLIB_GRIDALGORITHMCODE_INCLUDE_GUARD
#define _CODE_LIBS_TWODLIBLIB_GRIDALGORITHMCODE_INCLUDE_GUARD

#include <boost/filesystem.hpp>
#include <MPILib/include/utilities/Log.hpp>
#include <fstream>
#include <iostream>
#include "GridAlgorithm.hpp"
#include "Stat.hpp"
#include "TwoDLibException.hpp"
#include "display.hpp"
#include "GridReport.hpp"

namespace TwoDLib {

	GridAlgorithm::GridAlgorithm
	(
		const std::string& model_name,
		const std::string& transform_matrix,
		MPILib::Time h,
		double start_v,
		double start_w,
		MPILib::Time tau_refractive,
		const std::string&  rate_method
	):
	_model_name(model_name),
	_rate_method(rate_method),
	_rate(0.0),
	_t_cur(0.0),
	_root(CreateRootNode(model_name)),
	_vec_mesh(CreateMeshObject()),
	_vec_vec_rev(std::vector<std::vector<Redistribution> >{this->Mapping("Reversal")}),
	_vec_vec_res(std::vector<std::vector<Redistribution> >{this->Mapping("Reset")}),
	_vec_tau_refractive(std::vector<MPILib::Time>({tau_refractive})),
	_dt(_vec_mesh[0].TimeStep()),
	_sys(_vec_mesh,_vec_vec_rev,_vec_vec_res,_vec_tau_refractive),
	_n_evolve(0),
	_n_steps(0),
	_sysfunction(rate_method == "AvgV" ? &TwoDLib::Ode2DSystemGroup::AvgV : &TwoDLib::Ode2DSystemGroup::F),
	_start_v(start_v),
	_start_w(start_w),
	_transform_matrix(transform_matrix)
	{
		_mass_swap = vector<double>(_sys._vec_mass.size());

		vector<Coordinates> coords = _vec_mesh[0].findPointInMeshSlow(Point(start_v, start_w));

		_sys.Initialize(0,coords[0][0],coords[0][1]);

		for(int i=0; i< 10000; i++)
			_sys._individuals[i] = _sys.Map(0,coords[0][0],coords[0][1]);

		_sys.setIndividuals(_sys._individuals);

	}

	GridAlgorithm::GridAlgorithm(const GridAlgorithm& rhs):
	_model_name(rhs._model_name),
	_rate_method(rhs._rate_method),
	_rate(rhs._rate),
	_t_cur(rhs._t_cur),
	_vec_mesh(rhs._vec_mesh),
	_root(rhs._root),
	_vec_vec_rev(rhs._vec_vec_rev),
	_vec_vec_res(rhs._vec_vec_res),
	_dt(_vec_mesh[0].TimeStep()),
	_vec_tau_refractive(rhs._vec_tau_refractive),
	_sys(_vec_mesh,_vec_vec_rev,_vec_vec_res,_vec_tau_refractive),
	_n_evolve(0),
	_n_steps(0),
	_sysfunction(rhs._sysfunction),
	_start_v(rhs._start_v),
	_start_w(rhs._start_w),
	_transform_matrix(rhs._transform_matrix)
	{
		_mass_swap = vector<double>(_sys._vec_mass.size());

		vector<Coordinates> coords = _vec_mesh[0].findPointInMeshSlow(Point(_start_v, _start_w));

		_sys.Initialize(0,coords[0][0],coords[0][1]);

		for(int i=0; i< 10000; i++)
			_sys._individuals[i] = _sys.Map(0,coords[0][0],coords[0][1]);

		_sys.setIndividuals(_sys._individuals);

	}

	GridAlgorithm* GridAlgorithm::clone() const
	{
	  return new GridAlgorithm(*this);
	}

	void GridAlgorithm::assignNodeId( MPILib::NodeId nid ) {
		_node_id = nid;
	}

	void GridAlgorithm::configure(const MPILib::SimulationRunParameter& par_run)
	{
		_transformMatrix = TransitionMatrix(_transform_matrix);
		_csr_transform = new CSRMatrix(_transformMatrix, _sys, 0);

		Display::getInstance()->addOdeSystem(_node_id, &_sys);
		GridReport<CustomConnectionParameters>::getInstance()->registerObject(_node_id, this);

		_t_cur = par_run.getTBegin();
		_network_time_step     = par_run.getTStep();

		_sys.InitializeResetRefractive(_network_time_step);

		Quadrilateral q1 = _sys.MeshObjects()[0].Quad(1,0);
		Quadrilateral q2 = _sys.MeshObjects()[0].Quad(1,1);

		double cell_h_dist = std::fabs(q2.Centroid()[0] - q1.Centroid()[0]);
		double cell_v_dist = std::fabs(q2.Centroid()[1] - q1.Centroid()[1]);

		// one of these distances should be close to zero, so pick the other one
		// we do this because we don't know if this is a v- or h- efficacy

		double cell_width = std::max(cell_h_dist, cell_v_dist);

		setupMasterSolver(cell_width);
		// at this stage initialization must have taken place, either by default in (0,0),
		// or by the user calling Initialize if there is no strip 0

		double sum = _sys.P();
		if (sum == 0.)
			throw TwoDLib::TwoDLibException("No initialization of the mass array has taken place. Call Initialize before configure.");

	}

	void GridAlgorithm::setupMasterSolver(double cell_width){
		try {
			std::unique_ptr<MasterGrid> p_master(new MasterGrid(_sys,cell_width));
			_p_master = std::move(p_master);
		}
		// TODO: investigate the following
		// for some reason, the exception is usually not caught by the main program, which is why we write its message to cerr here.
		catch(TwoDLibException& e){
			std::cerr << e.what() << std::endl;
			throw e;
		}

	}

	std::vector<TwoDLib::Redistribution> GridAlgorithm::Mapping(const string& type)
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

	std::vector<Mesh> GridAlgorithm::CreateMeshObject(){
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
		std::vector<TwoDLib::Mesh> vec_mesh{ mesh };
		return vec_mesh;
	}

	pugi::xml_node GridAlgorithm::CreateRootNode(const string& model_name){

		// document
		pugi::xml_parse_result result = _doc.load_file(model_name.c_str());
		pugi::xml_node  root = _doc.first_child();

		if (result.status != pugi::status_ok)
		  throw TwoDLib::TwoDLibException("Can't open .model file.");
		return root;
	}

	void GridAlgorithm::evolveNodeState
	(
		const std::vector<MPILib::Rate>& nodeVector,
		const std::vector<CustomConnectionParameters>& weightVector,
		MPILib::Time time
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
			  LOG(MPILib::utilities::logWARNING) << "Mesh runs at a time step which is a multiple of the network time step. Is this intended?";
			else
			  ; // else is fine
		}
	    // mass rotation
	    for (MPILib::Index i = 0; i < _n_steps; i++){

				_sys.EvolveWithoutMeshUpdate();
#pragma omp parallel for
				 for(unsigned int id = 0; id < _mass_swap.size(); id++)
					  _mass_swap[id] = 0.;

			 	_csr_transform->MV(_mass_swap,_sys._vec_mass);
				// _csr_transform->MVIndividual(_sys._individuals);

				_sys._vec_mass = _mass_swap;
	    }

			// WARNING: originally reset goes afvirtual void applyMasterSolver();ter master but this way,
			// we can guarantee there's no mass above threshold when running
			// MVGrid in MasterGrid
			_sys.RedistributeProbability(_n_steps);
			// _sys.RedistributeIndividual(350, 150, 30, _n_steps);

			std::vector<std::vector<MPILib::Rate>> vec_rates;
			for(unsigned int i=0; i<_vec_vec_delay_queues.size(); i++){
				std::vector<MPILib::Rate> rates;
				for(unsigned int j=0; j<_vec_vec_delay_queues[i].size(); j++){
					rates.push_back(_vec_vec_delay_queues[i][j].getCurrentRate());
				}
				vec_rates.push_back(rates);
			}

			applyMasterSolver(vec_rates[0]);

			_t_cur += _n_steps*_dt;

 	    _rate = (_sys.*_sysfunction)()[0];

 	    _n_evolve++;
	}

	void GridAlgorithm::applyMasterSolver(std::vector<MPILib::Rate> rates) {
			_p_master->Apply(_n_steps*_dt,rates,_efficacy_map);
			// _p_master->ApplyIndividual(_n_steps*_dt,rates,_efficacy_map, _sys._individuals);
			// _p_master->ApplyIndividualPoisson(_n_steps*_dt,rates,_efficacy_map, _sys._individuals);
	}

	void GridAlgorithm::FillMap(const std::vector<CustomConnectionParameters>& vec_weights)
	{
		// this function will only be called once;
		_efficacy_map = std::vector<double>(vec_weights.size());

 		for(MPILib::Index i_weight = 0; i_weight < _efficacy_map.size(); i_weight++){
			_efficacy_map[i_weight] = std::stod(vec_weights[i_weight]._params.at("efficacy"));
		}

		_vec_vec_delay_queues = std::vector< std::vector<MPILib::DelayedConnectionQueue> >(0); // MeshAlgorithm really only uses the first array, i.e. the rates it receives in prepareEvole
 		_vec_vec_delay_queues.push_back( std::vector<MPILib::DelayedConnectionQueue>(vec_weights.size()));
		for(unsigned int q = 0; q < vec_weights.size(); q++){
			_vec_vec_delay_queues[0][q] = MPILib::DelayedConnectionQueue(_network_time_step, std::stod(vec_weights[q]._params.at("delay")));
		}
	}

	MPILib::AlgorithmGrid GridAlgorithm::getGrid(MPILib::NodeId id, bool b_state) const
	{
		// An empty grid will lead to crashes
		vector<double> array_interpretation {0.};
		vector<double> array_state {0.};

		return MPILib::AlgorithmGrid(array_state,array_interpretation);
	}

	void GridAlgorithm::reportDensity(MPILib::Time t) const
	{
		std::ostringstream ost;
		ost << _node_id  << "_" << t;
		ost << "_" << _sys.P();
		string fn("density_mesh_" + ost.str());

		std::string model_path = _model_name;
		boost::filesystem::path path(model_path);

		// MdK 27/01/2017. grid file is now created in the cwd of the program and
		// not in the directory where the mesh resides.
		const std::string dirname = path.filename().string() + "_mesh";

		if (! boost::filesystem::exists(dirname) ){
			boost::filesystem::create_directory(dirname);
		}
		std::ofstream ofst(dirname + "/" + fn);
		std::vector<std::ostream*> vec_str{&ofst};
		_sys.Dump(vec_str);
	}

	void GridAlgorithm::prepareEvolve
	(
		const std::vector<MPILib::Rate>& nodeVector,
		const std::vector<CustomConnectionParameters>& weightVector,
		const std::vector<MPILib::NodeType>& typeVector
	)
	{
		if (_efficacy_map.size() == 0){
			FillMap(weightVector);
		}

		// take into account the number of connections

		assert(nodeVector.size() == weightVector.size());
		for (MPILib::Index i = 0; i < nodeVector.size(); i++){
			double offset = 0.0;
			if (weightVector[i]._params.find("avgv_offset") != weightVector[i]._params.end())
				offset = std::stod(weightVector[i]._params.at("avgv_offset"));

			_vec_vec_delay_queues[0][i].updateQueue((offset + nodeVector[i])*std::stod(weightVector[i]._params.at("num_connections")));
		}

	}


}

#endif
