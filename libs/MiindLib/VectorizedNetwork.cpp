#include "VectorizedNetwork.hpp"
#include <boost/timer/timer.hpp>

using namespace MiindLib;

VectorizedNetwork::VectorizedNetwork(MPILib::Time time_step):
_num_nodes(0),
_n_steps(10),
_network_time_step(time_step)
{
}

void VectorizedNetwork::addGridNode(TwoDLib::Mesh mesh, TwoDLib::TransitionMatrix tmat, double start_v, double start_w,
  std::vector<TwoDLib::Redistribution> vec_rev, std::vector<TwoDLib::Redistribution> vec_res) {

  _num_nodes++;
  _node_id_to_group_mesh.insert(std::pair<MPILib::NodeId, MPILib::Index>(_num_nodes-1,_vec_mesh.size()));
  _group_mesh_to_node_id.insert(std::pair<MPILib::Index, MPILib::NodeId>(_vec_mesh.size(),_num_nodes-1));
  _grid_meshes.push_back(_vec_mesh.size());
  _vec_mesh.push_back(mesh);
  _vec_transforms.push_back(tmat);
  _start_vs.push_back(start_v);
  _start_ws.push_back(start_w);
  _vec_vec_rev.push_back(vec_rev);
  _vec_vec_res.push_back(vec_res);

}

void VectorizedNetwork::addMeshNode(TwoDLib::Mesh mesh,
  std::vector<TwoDLib::Redistribution> vec_rev, std::vector<TwoDLib::Redistribution> vec_res) {

  _num_nodes++;
  _node_id_to_group_mesh.insert(std::pair<MPILib::NodeId, MPILib::Index>(_num_nodes-1,_vec_mesh.size()));
  _group_mesh_to_node_id.insert(std::pair<MPILib::Index, MPILib::NodeId>(_vec_mesh.size(),_num_nodes-1));
  _mesh_meshes.push_back(_vec_mesh.size());
  _vec_mesh.push_back(mesh);
  _vec_vec_rev.push_back(vec_rev);
  _vec_vec_res.push_back(vec_res);

}

void VectorizedNetwork::addRateNode(function_pointer functor){
  _num_nodes++;
  _rate_functions.push_back(function_association(_num_nodes-1,functor));
}

void VectorizedNetwork::initOde2DSystem(){

  _group = new TwoDLib::Ode2DSystemGroup(_vec_mesh,_vec_vec_rev,_vec_vec_res);

	for( MPILib::Index i=0; i < _grid_meshes.size(); i++){
    vector<TwoDLib::Coordinates> coords = _vec_mesh[i].findPointInMeshSlow(TwoDLib::Point(_start_vs[i], _start_ws[i]));
    _group->Initialize(i,coords[0][0],coords[0][1]);

    //setup initial working index space for each grid mesh
    _current_index.push_back({_group->Map(i, coords[0][0],coords[0][1])-_group->Offsets()[i]});
    _current_indices_in_mesh.push_back(1);

    //create CSR Matrix for each transforms
    _csrs.push_back(TwoDLib::CSRMatrix(_vec_transforms[i], *(_group), i));
  }

  // All grids/meshes must have the same timestep
  TwoDLib::MasterParameter par(static_cast<MPILib::Number>(ceil(_vec_mesh[0].TimeStep()/_network_time_step)));
  _n_steps = std::max((int)par._N_steps,10);

  _group_adapter = new CudaTwoDLib::CudaOde2DSystemAdapter(*(_group));
}

void VectorizedNetwork::precalcWorkingIndexes(std::vector<inttype>& off1s, std::vector<inttype>& off2s){

  // calculate the largest noise spread
  inttype max_offset = 0;
  // find all cells which will be impacted by the master equation solver
  for(unsigned int o1=0; o1 < off1s.size(); o1++){
    max_offset = std::max((int)max_offset, (int)std::ceil((std::abs(off1s[o1]))));
  }
  for(unsigned int o2=0; o2 < off2s.size(); o2++){
    max_offset = std::max((int)max_offset, (int)std::ceil((std::abs(off2s[o2]))));
  }

  int noise_spread = max_offset*_n_steps;

  for (unsigned int m=0; m<_grid_meshes.size(); m++){

    unsigned int mesh_size = 0;
    for(unsigned int i=0; i<_vec_mesh[m].NrStrips();i++)
      for(unsigned int j=0; j<_vec_mesh[m].NrCellsInStrip(i);j++)
        mesh_size++;

    unsigned int mesh_offset = _group->Offsets()[m];

    unsigned int total_mass_size = _group->Mass().size();

    std::map<MPILib::Index, std::set<MPILib::Index>> map;

    for(unsigned int c=0; c<mesh_size; c++){
      std::set< MPILib::Index > c_spread;
      for(unsigned int t=0; t < _vec_transforms[m].Matrix()[c]._vec_to_line.size(); t++){
        unsigned int index = _group->Map(m,
          _vec_transforms[m].Matrix()[c]._vec_to_line[t]._to[0],
          _vec_transforms[m].Matrix()[c]._vec_to_line[t]._to[1])-mesh_offset;

        c_spread.insert(index);

          // expand the index to include all noise receiving cells
        for(int n=-noise_spread; n<noise_spread; n++){
          int new_ind = ((((int)c + (int)n)%(int)total_mass_size)+(int)total_mass_size)%(int)total_mass_size;
          c_spread.insert(new_ind);
        }
      }
      map.insert(std::pair<MPILib::Index,std::set<MPILib::Index>>(c,c_spread));
    } // end for each cell in mesh

    for(unsigned int c=0; c<_vec_vec_res[m].size(); c++){
      unsigned int c_full_index = _group->Map(m, _vec_vec_res[m][c]._from[0], _vec_vec_res[m][c]._from[1]);
      unsigned int reset_index = _group->Map(m, _vec_vec_res[m][c]._to[0], _vec_vec_res[m][c]._to[1]);

      for(int n=-noise_spread; n<noise_spread; n++){
        int new_ind = ((((int)reset_index + (int)n)%(int)total_mass_size)+(int)total_mass_size)%(int)total_mass_size;
        map[c_full_index-mesh_offset].insert(new_ind-mesh_offset);
      }
    } //end for each reset cell

    _grid_spread.push_back(map);
  } // for each mesh
}

void VectorizedNetwork::rectifyWorkingIndexes() {

  for(int m=0; m<_grid_meshes.size(); m++) {
    std::set<MPILib::Index> _new_current;
    std::vector<MPILib::Index>::iterator it;
    for(it = _current_index[m].begin(); it != _current_index[m].end(); it++){
      if (_group->Mass()[(*it)+_group->Offsets()[m]] > 0.0000001){
        _new_current.insert(*it);
        _new_current.insert(_grid_spread[m][*it].begin(), _grid_spread[m][*it].end());
      }
    }
    _current_index[m] = std::vector<MPILib::Index>(_new_current.begin(), _new_current.end());
    _current_indices_in_mesh[m] = _current_index[m].size();
  }
}

void VectorizedNetwork::reportNodeActivities(MPILib::Time sim_time){
  for (int i=0; i<_rate_nodes.size(); i++){
		std::ostringstream ost2;
		ost2 << "rate_" << i;
		std::ofstream ofst_rate(ost2.str(), std::ofstream::app);
		ofst_rate.precision(10);
		ofst_rate << sim_time << "\t" << _out_rates[i] << std::endl;
		ofst_rate.close();
	}
}

void VectorizedNetwork::addGridConnection(MPILib::NodeId in, MPILib::NodeId out, double efficacy, int n_conns){
  _grid_connections.push_back(NodeGridConnection(in,out,efficacy,n_conns));
}

void VectorizedNetwork::addMeshConnection(MPILib::NodeId in, MPILib::NodeId out, double efficacy, int n_conns, TwoDLib::TransitionMatrix *tmat){
  _mesh_connections.push_back(NodeMeshConnection(in,out,efficacy,n_conns,tmat));
}

void VectorizedNetwork::mainLoop(MPILib::Time t_begin, MPILib::Time t_end, MPILib::Time t_report, bool write_displays) {
	MPILib::Number n_iter = static_cast<MPILib::Number>(ceil((t_end - t_begin)/_network_time_step));
	MPILib::Number n_report = static_cast<MPILib::Number>(ceil((t_report - t_begin)/_network_time_step));

  for(unsigned int i=0; i<_display_nodes.size(); i++){
    TwoDLib::Display::getInstance()->addOdeSystem(i, _group, _node_id_to_group_mesh[i]);
  }

  const MPILib::Time h = 1./_n_steps*_vec_mesh[0].TimeStep();

  // Setup the OpenGL displays (if there are any required)
	TwoDLib::Display::getInstance()->animate(write_displays, _display_nodes, _network_time_step);

  // Generate calculated transition vectors for grid derivative
  std::vector<inttype> node_to_group_meshes;
  std::vector<fptype> stays;
  std::vector<fptype> goes;
  std::vector<inttype> off1s;
  std::vector<inttype> off2s;

  for (unsigned int i=0; i<_grid_connections.size(); i++){
    // for each connection, which of group's meshes is being affected
    node_to_group_meshes.push_back(_node_id_to_group_mesh[_grid_connections[i]._out]);
    // the input rate comes from the node indexed by connection _in
    TwoDLib::Mesh::GridCellTransition cell_transition =
      _group->MeshObjects()[_node_id_to_group_mesh[_grid_connections[i]._out]]
      .calculateCellTransition(_grid_connections[i]._efficacy);
    stays.push_back(cell_transition._stays);
    goes.push_back(cell_transition._goes);
    off1s.push_back(cell_transition._offset_1);
    off2s.push_back(cell_transition._offset_2);
  }

  precalcWorkingIndexes(off1s, off2s);

  for (unsigned int i=0; i<_mesh_connections.size(); i++){
    node_to_group_meshes.push_back(_node_id_to_group_mesh[_mesh_connections[i]._out]);
    _csrs.push_back(TwoDLib::CSRMatrix(*(_mesh_connections[i]._transition), *(_group), _node_id_to_group_mesh[_mesh_connections[i]._out]));
  }

  CudaTwoDLib::CSRAdapter _csr_adapter(*_group_adapter,_csrs,_grid_meshes.size(),_grid_connections.size()+_mesh_connections.size(),h);

  MPILib::utilities::ProgressBar *pb = new MPILib::utilities::ProgressBar(n_iter);
	MPILib::Time time = 0;
  boost::timer::auto_cpu_timer timer;
	for(MPILib::Index i_loop = 0; i_loop < n_iter; i_loop++){
		time = _network_time_step*i_loop;

    std::vector<fptype> rates;
    for (unsigned int i=0; i<_grid_connections.size(); i++){
      rates.push_back(_out_rates[_grid_connections[i]._in]*_grid_connections[i]._n_connections);
    }
    for (unsigned int i=0; i<_mesh_connections.size(); i++){
      rates.push_back(_out_rates[_mesh_connections[i]._in]*_mesh_connections[i]._n_connections);
    }

    rectifyWorkingIndexes();
    _group_adapter->transferGridIndex(_current_index);

		_group_adapter->Evolve(_mesh_meshes);

    _csr_adapter.ClearDerivative();
    // _csr_adapter.SingleTransformStep();
    _csr_adapter.SingleTransformStepIndexed(_current_indices_in_mesh);
    _csr_adapter.AddDerivativeFullIndexed(_current_indices_in_mesh);

    _group_adapter->ClearDerivative();
    _group_adapter->RedistributeProbabilityThreadedIndexed(_current_indices_in_mesh);
    _group_adapter->AddDerivativeFullIndexed(_current_indices_in_mesh);

    _group_adapter->MapFinishThreaded();

		for (MPILib::Index i_part = 0; i_part < _n_steps; i_part++ ){
			_csr_adapter.ClearDerivative();
      _csr_adapter.CalculateMeshGridDerivativeIndexed(node_to_group_meshes, rates, stays, goes, off1s, off2s, _current_indices_in_mesh);
			_csr_adapter.AddDerivative();
		}

    const std::vector<fptype>& group_rates = _group_adapter->F();
    for(unsigned int i=0; i<group_rates.size(); i++){
      _out_rates[_group_mesh_to_node_id[i]] = group_rates[i];
    }

    for( const auto& element: _rate_functions){
  		_out_rates[element.first] = element.second(time);
    }

    _group_adapter->updateGroupMass();
    TwoDLib::Display::getInstance()->updateDisplay(i_loop);
		reportNodeActivities(time);

    (*pb)++;

	}
}
