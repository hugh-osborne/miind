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

#ifndef CODE_MPILIB_MPINODE_HPP_
#define CODE_MPILIB_MPINODE_HPP_

#include <MPILib/include/MPINode.hpp>
#include <iostream>
#include <MPILib/include/BasicDefinitions.hpp>
#include <MPILib/include/utilities/Exception.hpp>
#include <MPILib/include/utilities/Log.hpp>
#include <MPILib/include/utilities/MPIProxy.hpp>
namespace MPILib {

template<class Weight, class NodeDistribution>
MPINode<Weight, NodeDistribution>::MPINode(
		const AlgorithmInterface<Weight>& algorithm,
		NodeType nodeType,
		NodeId nodeId,
		const NodeDistribution& nodeDistribution,
		const std::map<NodeId, MPINode<Weight, NodeDistribution>>& localNode,
		const std::string& name) :
		_pAlgorithm(algorithm.clone()), //
		_nodeType(nodeType), //
		_nodeId(nodeId), //
		_rLocalNodes(localNode), //
		_rNodeDistribution(nodeDistribution),
		_name(name) {
}

template<class Weight, class NodeDistribution>
MPINode<Weight, NodeDistribution>::~MPINode() {
}

template<class Weight, class NodeDistribution>
Time MPINode<Weight, NodeDistribution>::evolve(Time time) {

	// A Node will call its Algorithm to update its state up until time 'time'. It will
	// return the time maintained by the Algorithm. This time may be sligthly different due to
	// rounding errors. Network::evolve will use this time to check whether Algorithms keep synchronized
	// with the overall network within reasonable bounds.

	// MdK: 05/07/2017. removed a while loop. The algorithm is now entirely responsible
	// for evolution up until the required time. The MPINetwork::evolve method will
	// now check whether algorithms keep consistent time.
	//printf("PROC %i evolved node %i before copy.\n", utilities::MPIProxy().getRank(), _nodeId);
	std::vector<ActivityType> _pActivity(_precursorActivity);
	std::vector<Weight> _pWeights(_weights);
	std::vector<NodeType> _pTypes(_precursorTypes);
	if(_hasExternalPrecursor) {
		_pActivity.push_back(_externalPrecursorActivity);
		_pWeights.push_back(_externalPrecursorWeight);
		_pTypes.push_back(_externalPrecursorType);
	}

	++_number_iterations;

	_pAlgorithm->evolveNodeState(_pActivity, _pWeights, time, _pTypes);
	Time t_ret = _pAlgorithm->getCurrentTime();

	if (fabs(t_ret - time) > MPILib::ALGORITHM_NETWORK_DISCREPANCY ){
		throw MPILib::utilities::Exception("There is a discrepancy between Algorithm and Network time");
	}

	// update state
	this->setActivity(_pAlgorithm->getCurrentRate());

	sendOwnActivity();
	receiveData();

	return _pAlgorithm->getCurrentTime();
}

template<class Weight, class NodeDistribution>
void MPINode<Weight, NodeDistribution>::prepareEvolve() {

	std::vector<ActivityType> _pActivity(_precursorActivity);
	std::vector<Weight> _pWeights(_weights);
	std::vector<NodeType> _pTypes(_precursorTypes);
	if(_hasExternalPrecursor) {
		_pActivity.push_back(_externalPrecursorActivity);
		_pWeights.push_back(_externalPrecursorWeight);
		_pTypes.push_back(_externalPrecursorType);
	}

	_pAlgorithm->prepareEvolve(_pActivity, _pWeights, _pTypes);

}

template<class Weight, class NodeDistribution>
ActivityType MPINode<Weight, NodeDistribution>::getActivity(){
	return _activity;
}

template<class Weight, class NodeDistribution>
void MPINode<Weight, NodeDistribution>::setExternalPrecurserActivity(ActivityType activity){
	_externalPrecursorActivity = activity;
}

template<class Weight, class NodeDistribution>
void MPINode<Weight, NodeDistribution>::recvExternalPrecurserActivity(NodeId id, int tag){
	utilities::MPIProxy().irecv(id, tag, _externalPrecursorActivity);
}

template<class Weight, class NodeDistribution>
void MPINode<Weight, NodeDistribution>::configureSimulationRun(
		const SimulationRunParameter& simParam) {

	_maximum_iterations = simParam.getMaximumNumberIterations();

	_pAlgorithm->assignNodeId(_nodeId);
	_pAlgorithm->configure(simParam);

	// Add this line or other nodes will not get a proper input at the first simulation step!
	this->setActivity(_pAlgorithm->getCurrentRate());

	_pHandler = std::shared_ptr<report::handler::AbstractReportHandler>(
			simParam.getHandler().clone());

	_pHandler->initializeHandler(_nodeId);

}

template<class Weight, class NodeDistribution>
void MPINode<Weight, NodeDistribution>::setExternalPrecursor(const Weight& weight, NodeType nodeType) {
			_hasExternalPrecursor = true;
			_externalPrecursorWeight = weight;
			_externalPrecursorType = nodeType;
}

template<class Weight, class NodeDistribution>
void MPINode<Weight, NodeDistribution>::addPrecursor(NodeId nodeId,
		const Weight& weight, NodeType nodeType) {
	_precursors.push_back(nodeId);
	_precursorTypes.push_back(nodeType);
	_weights.push_back(weight);
	//make sure that _precursorStates is big enough to store the data
	_precursorActivity.resize(_precursors.size());
}

template<class Weight, class NodeDistribution>
void MPINode<Weight, NodeDistribution>::addSuccessor(NodeId nodeId) {
	_successors.push_back(nodeId);
}

template<class Weight, class NodeDistribution>
ActivityType MPINode<Weight, NodeDistribution>::getActivity() const {
	return _activity;
}

template<class Weight, class NodeDistribution>
NodeId MPINode<Weight, NodeDistribution>::getNodeId() const {
	return _nodeId;
}

template<class Weight, class NodeDistribution>
ActivityType MPINode<Weight, NodeDistribution>::getExternalPrecursorActivity() {
	return _externalPrecursorActivity;
}

template<class Weight, class NodeDistribution>
void MPINode<Weight, NodeDistribution>::setActivity(ActivityType activity) {
	_activity = activity;
}

template<class Weight, class NodeDistribution>
void MPINode<Weight, NodeDistribution>::waitAll() {
	utilities::MPIProxy().waitAll();
}


template<class Weight, class NodeDistribution>
void MPINode<Weight, NodeDistribution>::receiveData() {
	int i = 0;

	for (auto it = _precursors.begin(); it != _precursors.end(); it++, i++) {
		//do not send the data if the node is local!
		if (_rNodeDistribution.isLocalNode(*it)) {
			_precursorActivity[i] =
					_rLocalNodes.find(*it)->second.getActivity();

		} else {
			utilities::MPIProxy().irecv(_rNodeDistribution.getResponsibleProcessor(*it), *it,
					_precursorActivity[i]);
		}
	}
}

template<class Weight, class NodeDistribution>
void MPINode<Weight, NodeDistribution>::sendOwnActivity() {

	for (auto& it : _successors) {
		//do not send the data if the node is local!
		if (!_rNodeDistribution.isLocalNode(it)) {
			utilities::MPIProxy().isend(_rNodeDistribution.getResponsibleProcessor(it),
					_nodeId, _activity);
		}
	}
}

template<class Weight, class NodeDistribution>
void MPINode<Weight, NodeDistribution>::reportAll(
		report::ReportType type) const {
	std::vector<report::ReportValue> vec_values;

	if (type == report::RATE) {
		report::Report report(_pAlgorithm->getCurrentTime(),
				Rate(this->getActivity()), this->_nodeId,
				_pAlgorithm->getGrid(this->_nodeId,false), type, vec_values, _rLocalNodes.size());
		_pHandler->writeReport(report);

	} else if ( type == report::STATE) {
		// We don't want getGrid to be called if there is no requirement to write the state, as this is expensive for some algorithms

		report::Report report(_pAlgorithm->getCurrentTime(),
				Rate(this->getActivity()), this->_nodeId,
				_pAlgorithm->getGrid(this->_nodeId,_pHandler->isStateWriteMandatory()), type, vec_values, _rLocalNodes.size());
		_pHandler->writeReport(report);
	}
}

template<class Weight, class NodeDistribution>
void MPINode<Weight, NodeDistribution>::clearSimulation() {

	_pHandler->detachHandler(_nodeId);
}

template<class Weight, class NodeDistribution>
NodeType MPINode<Weight, NodeDistribution>::getNodeType() const {
	return _nodeType;
}

}
//end namespace MPILib

#endif /* CODE_MPILIB_MPINODE_HPP_ */
