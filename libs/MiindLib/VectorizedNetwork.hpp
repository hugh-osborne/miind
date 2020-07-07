#ifndef _CODE_MIINDLIB_VectorizedNetwork_INCLUDE_GUARD
#define _CODE_MIINDLIB_VectorizedNetwork_INCLUDE_GUARD

#include <CudaTwoDLib/CudaTwoDLib.hpp>
#include <MPILib/include/DelayedConnectionQueue.hpp>
#include <TwoDLib/display.hpp>
#include <TwoDLib/MasterParameter.hpp>

typedef CudaTwoDLib::fptype fptype;
typedef CudaTwoDLib::inttype inttype;
typedef MPILib::Rate (*function_pointer)(MPILib::Time);
typedef std::pair<MPILib::Index, function_pointer> function_association;
typedef std::vector<function_association> function_list;

namespace MiindLib {

class NodeMeshConnection {
public:
  bool _external;
  unsigned int _external_id;
  MPILib::NodeId _in;
  MPILib::NodeId _out;
  double _efficacy;
  double _delay;
  TwoDLib::TransitionMatrix *_transition;
  int _n_connections;

  NodeMeshConnection(MPILib::NodeId in, MPILib::NodeId out, double eff, int n_conns, double delay, TwoDLib::TransitionMatrix *trans):
  _external(false),_external_id(0),_in(in),_out(out),_efficacy(eff),_n_connections(n_conns),_delay(delay),_transition(trans){}

  NodeMeshConnection(MPILib::NodeId out, double eff, int n_conns, double delay, TwoDLib::TransitionMatrix *trans, MPILib::NodeId ext_id):
  _external(true),_external_id(ext_id),_out(out),_efficacy(eff),_n_connections(n_conns),_delay(delay),_transition(trans){}
};

class NodeGridConnection {
public:
  bool _external;
  unsigned int _external_id;
  MPILib::NodeId _in;
  MPILib::NodeId _out;
  std::map<std::string, std::string> _params;

  NodeGridConnection(MPILib::NodeId in, MPILib::NodeId out, std::map<std::string, std::string> params):
  _external(false),_external_id(0),_in(in),_out(out),_params(params){}

  NodeGridConnection(MPILib::NodeId out, std::map<std::string, std::string> params, MPILib::NodeId ext_id):
  _external(true),_external_id(ext_id),_out(out),_params(params){}
};

class NodeMeshCustomConnection {
public:
  bool _external;
  unsigned int _external_id;
  MPILib::NodeId _in;
  MPILib::NodeId _out;
  TwoDLib::TransitionMatrix *_transition;
  std::map<std::string, std::string> _params;

  NodeMeshCustomConnection(MPILib::NodeId in, MPILib::NodeId out, std::map<std::string, std::string> params, TwoDLib::TransitionMatrix *trans):
  _external(false),_external_id(0),_in(in),_out(out),_params(params),_transition(trans){}

  NodeMeshCustomConnection(MPILib::NodeId out, std::map<std::string, std::string> params, TwoDLib::TransitionMatrix *trans, MPILib::NodeId ext_id):
  _external(true),_external_id(ext_id),_out(out),_params(params),_transition(trans){}
};

class VectorizedNetwork {
public:
  VectorizedNetwork(MPILib::Time time_step);

  void initOde2DSystem(unsigned int min_solve_steps=10);

  void setRateNodes(std::vector<MPILib::NodeId> ids, std::vector<MPILib::Time> intervals){
    _rate_nodes = ids;
    _rate_intervals = intervals;
  }
  void setDisplayNodes(std::vector<MPILib::NodeId> ids){
    _display_nodes = ids;
  }
  void setDensityNodes(std::vector<MPILib::NodeId> ids, std::vector<MPILib::Time> start_times,
    std::vector<MPILib::Time> end_times, std::vector<MPILib::Time> intervals){
    _density_nodes = ids;
    _density_start_times = start_times;
    _density_end_times = end_times;
    _density_intervals = intervals;
  }

  void addGridNode(TwoDLib::Mesh mesh, TwoDLib::TransitionMatrix tmat, double start_v, double start_w,
    std::vector<TwoDLib::Redistribution> vec_rev, std::vector<TwoDLib::Redistribution> vec_res, double tau_refractive);

  void addMeshNode(TwoDLib::Mesh mesh, std::vector<TwoDLib::Redistribution> vec_rev, std::vector<TwoDLib::Redistribution> vec_res, double tau_refractive);

  void addRateNode(function_pointer functor);

  void addExternalMonitor(MPILib::NodeId node);

  void addGridConnection(MPILib::NodeId in, MPILib::NodeId out, std::map<std::string, std::string> params);

  void addGridConnection(MPILib::NodeId out,std::map<std::string, std::string> params, MPILib::NodeId ext_id);

  void addMeshConnection(MPILib::NodeId in, MPILib::NodeId out, double efficacy, int n_conns, double delay, TwoDLib::TransitionMatrix *tmat);

  void addMeshConnection(MPILib::NodeId out, double efficacy, int n_conns, double delay, TwoDLib::TransitionMatrix *tmat, MPILib::NodeId ext_id);

  void addMeshCustomConnection(MPILib::NodeId in, MPILib::NodeId out, std::map<std::string, std::string> params, TwoDLib::TransitionMatrix *tmat);

  void addMeshCustomConnection(MPILib::NodeId out,std::map<std::string, std::string> params,  TwoDLib::TransitionMatrix *tmat, MPILib::NodeId ext_id);


  void reportNodeActivities(MPILib::Time sim_time);
  void reportNodeDensities(MPILib::Time sim_time);
  void mainLoop(MPILib::Time t_begin, MPILib::Time t_end, MPILib::Time t_report, bool write_displays);

  void setupLoop(bool write_displays);
  std::vector<double> singleStep(std::vector<double>, unsigned int i_loop);
  void endLoops();

protected:

  std::vector<TwoDLib::TransitionMatrix> _vec_transforms;
  std::vector<double> _start_vs;
  std::vector<double> _start_ws;

  std::vector<TwoDLib::Mesh> _vec_mesh;
  std::vector< std::vector<TwoDLib::Redistribution> > _vec_vec_rev;
  std::vector< std::vector<TwoDLib::Redistribution> > _vec_vec_res;
  std::vector<MPILib::Time> _vec_tau_refractive;

  std::vector<inttype> _grid_node_ids;
  std::vector<TwoDLib::Mesh> _grid_vec_mesh;
  std::vector< std::vector<TwoDLib::Redistribution> > _grid_vec_vec_rev;
  std::vector< std::vector<TwoDLib::Redistribution> > _grid_vec_vec_res;
  std::vector<MPILib::Time> _grid_vec_tau_refractive;

  std::vector<inttype> _mesh_node_ids;
  std::vector<TwoDLib::Mesh> _mesh_vec_mesh;
  std::vector< std::vector<TwoDLib::Redistribution> > _mesh_vec_vec_rev;
  std::vector< std::vector<TwoDLib::Redistribution> > _mesh_vec_vec_res;
  std::vector<MPILib::Time> _mesh_vec_tau_refractive;

  std::vector<inttype> _rate_func_node_ids;

  std::vector<inttype> _grid_transform_indexes;
  std::vector<inttype> _mesh_transform_indexes;

  std::vector<inttype> _grid_meshes;
  std::vector<inttype> _mesh_meshes;

  TwoDLib::Ode2DSystemGroup *_group;

  CudaTwoDLib::CudaOde2DSystemAdapter *_group_adapter;
  CudaTwoDLib::CSRAdapter *_csr_adapter;

  std::vector<TwoDLib::CSRMatrix> _csrs;

  unsigned int _n_steps;
  unsigned int _master_steps;

  unsigned int _num_nodes;
  MPILib::Time _network_time_step;

  std::vector<MPILib::NodeId> _display_nodes;
  std::vector<MPILib::NodeId> _rate_nodes;
  std::vector<MPILib::Time> _rate_intervals;
  std::vector<MPILib::NodeId> _density_nodes;
  std::vector<MPILib::Time> _density_start_times;
  std::vector<MPILib::Time> _density_end_times;
  std::vector<MPILib::Time> _density_intervals;

  std::vector<inttype> _connection_out_group_mesh;
  std::vector<MPILib::DelayedConnectionQueue> _connection_queue;
  std::map<MPILib::NodeId, std::vector<MPILib::NodeId>> _node_to_connection_queue;
  std::map<MPILib::NodeId, std::vector<MPILib::NodeId>> _external_to_connection_queue;
  std::vector<fptype> _effs;

  std::map<MPILib::NodeId, fptype> _current_node_rates;

  std::vector<inttype> _monitored_nodes;

  std::vector<NodeMeshConnection> _mesh_connections;
  std::vector<NodeMeshCustomConnection> _mesh_custom_connections;
  std::vector<NodeGridConnection> _grid_connections;

  function_list _rate_functions;

  std::map<MPILib::NodeId, MPILib::Index> _node_id_to_group_mesh;
  std::map<MPILib::Index, MPILib::NodeId> _group_mesh_to_node_id;

};

}

#endif
