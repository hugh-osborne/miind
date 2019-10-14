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

#include <iostream>
#include <cuda_runtime.h>
#include <cstdio>
#include <cmath>
#include "CudaEuler.cuh"
#include "CSRAdapter.cuh"

using namespace CudaTwoDLib;

const fptype TOLERANCE = 1e-9;


#define checkCudaErrors(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
  if (code != cudaSuccess) {
    fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    if (abort) exit(code);
  }
}

void CSRAdapter::FillMatrixMaps(const std::vector<TwoDLib::CSRMatrix>& vecmat)
{
   for(inttype m = 0; m < vecmat.size(); m++)
   {
       _nval[m] = vecmat[m].Val().size();
       checkCudaErrors(cudaMalloc((fptype**)&_val[m],_nval[m]*sizeof(fptype)));
       // dont't depend on Val() being of fptype
       std::vector<fptype> vecval;
       for (fptype val: vecmat[m].Val())
           vecval.push_back(val);
       checkCudaErrors(cudaMemcpy(_val[m],&vecval[0],sizeof(fptype)*_nval[m],cudaMemcpyHostToDevice));

       _nia[m] = vecmat[m].Ia().size();
       checkCudaErrors(cudaMalloc((inttype**)&_ia[m],_nia[m]*sizeof(inttype)));
       std::vector<inttype> vecia;
       for(inttype ia: vecmat[m].Ia())
           vecia.push_back(ia);
       checkCudaErrors(cudaMemcpy(_ia[m],&vecia[0],sizeof(inttype)*_nia[m],cudaMemcpyHostToDevice));


       _nja[m] = vecmat[m].Ja().size();
       checkCudaErrors(cudaMalloc((inttype**)&_ja[m],_nja[m]*sizeof(inttype)));
       std::vector<inttype> vecja;
       for(inttype ja: vecmat[m].Ja())
           vecja.push_back(ja);
       checkCudaErrors(cudaMemcpy(_ja[m],&vecja[0],sizeof(inttype)*_nja[m],cudaMemcpyHostToDevice));
   }
}

void CSRAdapter::InitializeStaticGridEfficacies(const std::vector<inttype>& vecindex,const std::vector<fptype>& efficacy) {
  _nr_grid_connections = efficacy.size();
  for(inttype m = 0; m < efficacy.size(); m++)
  {
    checkCudaErrors(cudaMalloc((fptype**)&_goes[m],_nr_rows[vecindex[m]]*sizeof(fptype)));
    checkCudaErrors(cudaMalloc((fptype**)&_stays[m],_nr_rows[vecindex[m]]*sizeof(fptype)));
    checkCudaErrors(cudaMalloc((inttype**)&_offset1s[m],_nr_rows[vecindex[m]]*sizeof(inttype)));
    checkCudaErrors(cudaMalloc((inttype**)&_offset2s[m],_nr_rows[vecindex[m]]*sizeof(inttype)));

    inttype numBlocks = (_nr_rows[vecindex[m]] + _blockSize - 1)/_blockSize;

    CudaCalculateGridEfficacies<<<numBlocks,_blockSize>>>(_nr_rows[vecindex[m]],
      efficacy[m], _cell_widths[vecindex[m]],
      _stays[m], _goes[m], _offset1s[m], _offset2s[m]);
  }
}

void CSRAdapter::InitializeStaticGridConductanceEfficacies(const std::vector<inttype>& vecindex,
  const std::vector<fptype>& efficacy, const std::vector<fptype>& rest_vs) {
    _nr_grid_connections = efficacy.size();

    checkCudaErrors(cudaMalloc((fptype**)&_cell_vs,_group.getGroup().Vs().size()*sizeof(fptype)));

    std::vector<fptype> vecval;
    for (double val: _group.getGroup().Vs())
        vecval.push_back((fptype)val);

    checkCudaErrors(cudaMemcpy(_cell_vs,&vecval[0],_group.getGroup().Vs().size()*sizeof(fptype),cudaMemcpyHostToDevice));

    for(inttype m = 0; m < efficacy.size(); m++)
    {
      checkCudaErrors(cudaMalloc((fptype**)&_goes[m],_nr_rows[vecindex[m]]*sizeof(fptype)));
      checkCudaErrors(cudaMalloc((fptype**)&_stays[m],_nr_rows[vecindex[m]]*sizeof(fptype)));
      checkCudaErrors(cudaMalloc((inttype**)&_offset1s[m],_nr_rows[vecindex[m]]*sizeof(inttype)));
      checkCudaErrors(cudaMalloc((inttype**)&_offset2s[m],_nr_rows[vecindex[m]]*sizeof(inttype)));

      inttype numBlocks = (_nr_rows[vecindex[m]] + _blockSize - 1)/_blockSize;

      CudaCalculateGridEfficaciesWithConductance<<<numBlocks,_blockSize>>>(_nr_rows[vecindex[m]],
        efficacy[m], _cell_widths[vecindex[m]], _cell_vs, rest_vs[m],
        _stays[m], _goes[m], _offset1s[m], _offset2s[m],_offsets[vecindex[m]]);
    }
  }


void CSRAdapter::DeleteMatrixMaps()
{
    for(inttype m = 0; m < _nr_m; m++)
    {
        cudaFree(_val[m]);
        cudaFree(_ia[m]);
        cudaFree(_ja[m]);
    }
}

inttype CSRAdapter::NumberIterations(const CudaOde2DSystemAdapter& group, fptype euler_timestep) const
{
    fptype tstep = group._group.MeshObjects()[0].TimeStep();
    for ( const auto& mesh: group._group.MeshObjects() )
        if (fabs(tstep - mesh.TimeStep()) > TOLERANCE){
           std::cerr << "Not all meshes in this group have the same time step. " <<  tstep << " " << mesh.TimeStep() << " " << tstep - mesh.TimeStep()  << std::endl;
           exit(0);
        }
    inttype  n_steps = static_cast<inttype>(std::round(tstep/euler_timestep));

    return n_steps;
}

void CSRAdapter::InspectMass(inttype i)
{
    std::vector<fptype> hostvec(_group._n);
    checkCudaErrors(cudaMemcpy(&hostvec[0],_group._mass,sizeof(fptype)*_group._n,cudaMemcpyDeviceToHost));
}

CSRAdapter::CSRAdapter(CudaOde2DSystemAdapter& group, const std::vector<TwoDLib::CSRMatrix>& vecmat,
 inttype nr_connections, fptype euler_timestep,
 const std::vector<inttype>& vecmat_indexes,const std::vector<inttype>& grid_transforms):
_group(group),
_euler_timestep(euler_timestep),
_nr_iterations(NumberIterations(group,euler_timestep)),
_nr_m(vecmat.size()),
_nr_streams(vecmat.size()),
_vecmats(vecmat_indexes),
_grid_transforms(grid_transforms),
_nval(std::vector<inttype>(vecmat.size())),
_val(std::vector<fptype*>(vecmat.size())),
_nia(std::vector<inttype>(vecmat.size())),
_ia(std::vector<inttype*>(vecmat.size())),
_nja(std::vector<inttype>(vecmat.size())),
_ja(std::vector<inttype*>(vecmat.size())),
_offsets(this->Offsets(vecmat)),
_nr_rows(this->NrRows(vecmat)),
_cell_widths(this->CellWidths(vecmat)),
_goes(std::vector<fptype*>(nr_connections)),
_stays(std::vector<fptype*>(nr_connections)),
_offset1s(std::vector<int*>(nr_connections)),
_offset2s(std::vector<int*>(nr_connections)),
_blockSize(256),
_numBlocks( (_group._n + _blockSize - 1) / _blockSize)
{
    this->FillMatrixMaps(vecmat);
    this->FillDerivative();
    this->CreateStreams();
}

CSRAdapter::CSRAdapter(CudaOde2DSystemAdapter& group, const std::vector<TwoDLib::CSRMatrix>& vecmat, fptype euler_timestep):
CSRAdapter(group,vecmat,vecmat.size(),euler_timestep,
std::vector<inttype>(),std::vector<inttype>()){
  for(unsigned int i=0; i<vecmat.size(); i++)
   _vecmats.push_back(i);
}

CSRAdapter::~CSRAdapter()
{
    this->DeleteMatrixMaps();
    this->DeleteDerivative();
    this->DeleteStreams();
}

void CSRAdapter::CreateStreams()
{
    _streams = (cudaStream_t *)malloc(_nr_streams*sizeof(cudaStream_t));
    for(int i = 0; i < _nr_streams; i++)
       cudaStreamCreate(&_streams[i]);
}

void CSRAdapter::DeleteStreams()
{
   free(_streams);
}

void CSRAdapter::FillDerivative()
{
    checkCudaErrors(cudaMalloc((fptype**)&_dydt,_group._n*sizeof(fptype)));
}

void CSRAdapter::DeleteDerivative()
{
    cudaFree(_dydt);
}

void CSRAdapter::ClearDerivative()
{
  inttype n=_group._n;
  CudaClearDerivative<<<_numBlocks,_blockSize>>>(n,_dydt,_group._mass);
}

std::vector<inttype> CSRAdapter::NrRows(const std::vector<TwoDLib::CSRMatrix>& vecmat) const
{
	std::vector<inttype> vecret;
	for (inttype m = 0; m < vecmat.size(); m++)
		vecret.push_back(vecmat[m].NrRows());
	return vecret;
}

std::vector<fptype> CSRAdapter::CellWidths(const std::vector<TwoDLib::CSRMatrix>& vecmat) const
{
	std::vector<fptype> vecret;
	for (inttype m = 0; m < _grid_transforms.size(); m++){
    vecret.push_back(_group.getGroup().MeshObjects()[_grid_transforms[m]].getCellWidth());
  }
	return vecret;
}

std::vector<inttype> CSRAdapter::Offsets(const std::vector<TwoDLib::CSRMatrix>& vecmat) const
{
	std::vector<inttype> vecret;
	for (inttype m = 0; m < vecmat.size(); m++)
		vecret.push_back(vecmat[m].Offset());
	return vecret;
}

void CSRAdapter::CalculateDerivative(const std::vector<fptype>& vecrates)
{
    for(inttype m : _vecmats)
    {
        // be careful to use this block size
        inttype numBlocks = (_nr_rows[m] + _blockSize - 1)/_blockSize;
        CudaCalculateDerivative<<<numBlocks,_blockSize>>>(_nr_rows[m],vecrates[m],_dydt,_group._mass,_val[m],_ia[m],_ja[m],_group._map,_offsets[m]);
    }

}

void CSRAdapter::CalculateGridDerivative(const std::vector<inttype>& vecindex, const std::vector<fptype>& vecrates, const std::vector<fptype>& vecstays, const std::vector<fptype>& vecgoes, const std::vector<int>& vecoff1s, const std::vector<int>& vecoff2s)
{
    for(inttype m = 0; m < vecindex.size(); m++)
    {
        // be careful to use this block size
        inttype numBlocks = (_nr_rows[vecindex[m]] + _blockSize - 1)/_blockSize;
        CudaCalculateGridDerivative<<<numBlocks,_blockSize,0,_streams[vecindex[m]]>>>(_nr_rows[vecindex[m]],vecrates[m],vecstays[m],vecgoes[m],vecoff1s[m],vecoff2s[m],_dydt,_group._mass,_offsets[m]);
    }

    cudaDeviceSynchronize();
}

void CSRAdapter::CalculateMeshGridDerivative(const std::vector<inttype>& vecindex,
  const std::vector<fptype>& vecrates, const std::vector<fptype>& vecstays,
  const std::vector<fptype>& vecgoes, const std::vector<int>& vecoff1s,
  const std::vector<int>& vecoff2s)
{

  for(inttype m = 0; m < vecstays.size(); m++)
  {
    // be careful to use this blo							void CalculateMeshGridDerivativeForward(const std::vector<inttype>& vecindex,ck size
    inttype numBlocks = (_nr_rows[vecindex[m]] + _blockSize - 1)/_blockSize;
    CudaCalculateGridDerivative<<<numBlocks,_blockSize,0,_streams[vecindex[m]]>>>(_nr_rows[vecindex[m]],vecrates[m],vecstays[m],vecgoes[m],vecoff1s[m],vecoff2s[m],_dydt,_group._mass,_offsets[vecindex[m]]);
  }

  for(int n=vecstays.size(); n<vecrates.size(); n++)
  {
    inttype mat_index = _grid_transforms.size() + (n - vecstays.size());
    // be careful to use this block size
    inttype numBlocks = (_nr_rows[mat_index] + _blockSize - 1)/_blockSize;
    CudaCalculateDerivative<<<numBlocks,_blockSize,0,_streams[vecindex[n]]>>>(_nr_rows[mat_index],vecrates[n],_dydt,_group._mass,_val[mat_index],_ia[mat_index],_ja[mat_index],_group._map,_offsets[mat_index]);
  }

  cudaDeviceSynchronize();

}

void CSRAdapter::CalculateMeshGridDerivativeWithEfficacy(const std::vector<inttype>& vecindex,
  const std::vector<fptype>& vecrates)
{
  for(inttype m = 0; m < _nr_grid_connections; m++)
  {
    // be careful to use this block size
    inttype numBlocks = (_nr_rows[vecindex[m]] + _blockSize - 1)/_blockSize;

    CudaCalculateGridDerivativeWithEfficacy<<<numBlocks,_blockSize,0,_streams[vecindex[m]]>>>(_nr_rows[vecindex[m]],
      vecrates[m],_stays[m], _goes[m], _offset1s[m], _offset2s[m],_dydt,_group._mass,_offsets[vecindex[m]]);
  }

  for(int n=_nr_grid_connections; n<vecrates.size(); n++)
  {
    inttype mat_index = _grid_transforms.size() + (n - _nr_grid_connections);
    // be careful to use this block size
    inttype numBlocks = (_nr_rows[mat_index] + _blockSize - 1)/_blockSize;
    CudaCalculateDerivative<<<numBlocks,_blockSize,0,_streams[vecindex[n]]>>>(_nr_rows[mat_index],vecrates[n],_dydt,_group._mass,_val[mat_index],_ia[mat_index],_ja[mat_index],_group._map,_offsets[mat_index]);
  }

  cudaDeviceSynchronize();

}

void CSRAdapter::SingleTransformStep()
{
  for(inttype m : _grid_transforms)
  {
      // be careful to use this block size
      inttype numBlocks = (_nr_rows[m] + _blockSize - 1)/_blockSize;
      CudaSingleTransformStep<<<numBlocks,_blockSize,0,_streams[m]>>>(_nr_rows[m],_dydt,_group._mass,_val[m],_ia[m],_ja[m],_group._map,_offsets[m]);
  }
}

void CSRAdapter::AddDerivative()
{
  EulerStep<<<_numBlocks,_blockSize>>>(_group._n,_dydt,_group._mass,_euler_timestep);
}

void CSRAdapter::AddDerivativeFull()
{
  EulerStep<<<_numBlocks,_blockSize>>>(_group._n,_dydt,_group._mass, 1.0);
}
