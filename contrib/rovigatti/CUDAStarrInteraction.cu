/*
 * CUDAStarrInteraction.cu
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#include "CUDAStarrInteraction.h"

#include "../Lists/CUDASimpleVerletList.h"
#include "../Lists/CUDANoList.h"
#include "../Particles/CustomParticle.h"

/* CUDA constants */
__constant__ int MD_N[1];
__constant__ int MD_N_hubs[1];
__constant__ float MD_box_side[1];

__constant__ float MD_LJ_sigma[3];
__constant__ float MD_LJ_sqr_sigma[3];
__constant__ float MD_LJ_sqr_rcut[3];
__constant__ float MD_LJ_E_cut[3];
__constant__ float MD_der_LJ_E_cut[3];
__constant__ float MD_fene_K[1];
__constant__ float MD_fene_sqr_r0[1];
__constant__ float MD_lin_k[1];
__constant__ float MD_sqr_rcut[1];

#include "../cuda_utils/CUDA_lr_common.cuh"

template <typename number, typename number4>
__device__ number4 minimum_image(const number4 &r_i, const number4 &r_j) {
	number dx = r_j.x - r_i.x;
	number dy = r_j.y - r_i.y;
	number dz = r_j.z - r_i.z;

	dx -= floorf(dx/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];
	dy -= floorf(dy/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];
	dz -= floorf(dz/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];

	return make_number4<number, number4>(dx, dy, dz, (number) 0.f);
}

template <typename number, typename number4>
__device__ number quad_minimum_image_dist(const number4 &r_i, const number4 &r_j) {
	number dx = r_j.x - r_i.x;
	number dy = r_j.y - r_i.y;
	number dz = r_j.z - r_i.z;

	dx -= floorf(dx/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];
	dy -= floorf(dy/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];
	dz -= floorf(dz/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];

	return dx*dx + dy*dy + dz*dz;
}

template <typename number, typename number4>
__device__ void _two_body(number4 &r, int int_type, int can_form_hb, number4 &F) {
	if(int_type == 2 && !can_form_hb) int_type = 1;

	number sqr_r = CUDA_DOT(r, r);

	number part = powf(MD_LJ_sqr_sigma[int_type]/sqr_r, 3.f);
	number force_mod = 24.f * part * (2.f*part - 1.f) / sqr_r;

	if(sqr_r > MD_LJ_sqr_rcut[int_type]) force_mod = (number) 0.f;

	F.x -= r.x * force_mod;
	F.y -= r.y * force_mod;
	F.z -= r.z * force_mod;
}

template<typename number, typename number4>
__device__ void _fene(number4 &r, number4 &F) {
	number sqr_r = CUDA_DOT(r, r);
	// this number is the module of the force over r, so we don't have to divide the distance
	// vector by its module
	number force_mod = -MD_fene_K[0] * MD_fene_sqr_r0[0] / (MD_fene_sqr_r0[0] - sqr_r);
	F.x -= r.x * force_mod;
	F.y -= r.y * force_mod;
	F.z -= r.z * force_mod;
}

template <typename number, typename number4>
__device__ void _particle_particle_bonded_interaction(number4 &ppos, number4 &qpos, number4 &F) {
	int ptype = get_particle_type<number, number4>(ppos);
	int qtype = get_particle_type<number, number4>(qpos);
	int int_type = ptype + qtype;

//	int pbtype = get_particle_btype<number, number4>(ppos);
//	int qbtype = get_particle_btype<number, number4>(qpos);
//	int int_btype = pbtype + qbtype;

	number4 r = qpos - ppos;
	_two_body<number, number4>(r, int_type, false, F);
	_fene<number, number4>(r, F);
}

template <typename number, typename number4>
__device__ void _particle_particle_interaction(number4 &ppos, number4 &qpos, number4 &F) {
	int ptype = get_particle_type<number, number4>(ppos);
	int qtype = get_particle_type<number, number4>(qpos);
	int int_type = ptype + qtype;

	int pbtype = get_particle_btype<number, number4>(ppos);
	int qbtype = get_particle_btype<number, number4>(qpos);
	int int_btype = pbtype + qbtype;
	int can_form_hb = (int_btype == 3);

	number4 r = minimum_image<number, number4>(ppos, qpos);
	_two_body<number, number4>(r, int_type, can_form_hb, F);
}

template<typename number, typename number4>
__device__ void _three_body(number4 &ppos, LR_bonds &bs, number4 &F, number4 *poss, number4 *forces) {
	if(bs.n3 == P_INVALID || bs.n5 == P_INVALID) return;

	number4 n3_pos = poss[bs.n3];
	number4 n5_pos = poss[bs.n5];

	number4 dist_pn3 = n3_pos - ppos;
	number4 dist_pn5 = ppos - n5_pos;

	number sqr_dist_pn3 = CUDA_DOT(dist_pn3, dist_pn3);
	number sqr_dist_pn5 = CUDA_DOT(dist_pn5, dist_pn5);
	number i_pn3_pn5 = 1.f / sqrtf(sqr_dist_pn3*sqr_dist_pn5);
	number cost = CUDA_DOT(dist_pn3, dist_pn5) * i_pn3_pn5;

	number cost_n3 = cost / sqr_dist_pn3;
	number cost_n5 = cost / sqr_dist_pn5;
	number force_mod_n3 = i_pn3_pn5 + cost_n3;
	number force_mod_n5 = i_pn3_pn5 + cost_n5;

	F += dist_pn3*(force_mod_n3*MD_lin_k[0]) - dist_pn5*(force_mod_n5*MD_lin_k[0]);
	forces[bs.n3] -= dist_pn3*(cost_n3*MD_lin_k[0]) - dist_pn5*(i_pn3_pn5*MD_lin_k[0]);
	forces[bs.n5] -= dist_pn3*(i_pn3_pn5*MD_lin_k[0]) - dist_pn5*(cost_n5*MD_lin_k[0]);
}

// forces + second step without lists
template <typename number, typename number4>
__global__ void Starr_forces(number4 *poss, number4 *forces, LR_bonds *bonds) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	LR_bonds bs = bonds[IND];
	number4 ppos = poss[IND];

	if(bs.n3 != P_INVALID) {
		number4 qpos = poss[bs.n3];
		_particle_particle_bonded_interaction<number, number4>(ppos, qpos, F);
	}

	if(bs.n5 != P_INVALID) {
		number4 qpos = poss[bs.n5];
		_particle_particle_bonded_interaction<number, number4>(ppos, qpos, F);
	}

	_three_body<number, number4>(ppos, bs, F, poss, forces);

	for(int j = 0; j < MD_N[0]; j++) {
		if(j != IND && bs.n3 != j && bs.n5 != j) {
			number4 qpos = poss[j];
			_particle_particle_interaction<number, number4>(ppos, qpos, F);
		}
	}

	// the real energy per particle is half of the one computed (because we count each interaction twice)
	F.w *= (number) 0.5f;
	forces[IND] = F;
}

// forces + second step with verlet lists
template <typename number, typename number4>
__global__ void Starr_forces(number4 *poss, number4 *forces, int *matrix_neighs, int *number_neighs, LR_bonds *bonds) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	number4 ppos = poss[IND];
	LR_bonds bs = bonds[IND];

	if(bs.n3 != P_INVALID) {
		number4 qpos = poss[bs.n3];
		_particle_particle_bonded_interaction<number, number4>(ppos, qpos, F);
	}
	if(bs.n5 != P_INVALID) {
		number4 qpos = poss[bs.n5];
		_particle_particle_bonded_interaction<number, number4>(ppos, qpos, F);
	}

	_three_body<number, number4>(ppos, bs, F, poss, forces);

	const int num_neighs = number_neighs[IND];
	for(int j = 0; j < num_neighs; j++) {
		const int k_index = matrix_neighs[j*MD_N[0] + IND];

		number4 qpos = poss[k_index];
		_particle_particle_interaction<number, number4>(ppos, qpos, F);
	}

	// the real energy per particle is half the one computed (because we count each interaction twice)
	F.w *= (number) 0.5f;

	forces[IND] = F;
}

template <typename number, typename number4>
__global__ void tetra_hub_forces(number4 *poss, number4 *forces, int *hubs, tetra_hub_bonds *bonds) {
	if(IND >= MD_N_hubs[0]) return;

	int idx_hub = hubs[IND];
	tetra_hub_bonds hub_bonds = bonds[IND];
	number4 poss_hub = poss[idx_hub];
	number4 F = forces[idx_hub];

	for(int an = 0; an < HUB_SIZE; an++) {
		int bonded_neigh = hub_bonds.n[an];
		// since bonded neighbours of anchors are in the anchor's neighbouring list, the LJ interaction between
		// the two, from the point of view of the anchor, has been already computed and hence the anchor-particle
		// interaction reduces to just the fene
		number4 r = poss[bonded_neigh] - poss_hub;
		_fene<number, number4>(r, F);
	}

	forces[idx_hub] = F;
}

template<typename number, typename number4>
CUDAStarrInteraction<number, number4>::CUDAStarrInteraction() : _h_tetra_hubs(NULL), _h_tetra_hub_neighs(NULL) {
	_N_hubs = -1;
	_d_tetra_hubs = NULL;
	_d_tetra_hub_neighs = NULL;
}

template<typename number, typename number4>
CUDAStarrInteraction<number, number4>::~CUDAStarrInteraction() {
	if(_h_tetra_hubs != NULL) {
		delete[] _h_tetra_hubs;
		delete[] _h_tetra_hub_neighs;

		CUDA_SAFE_CALL( cudaFree(_d_tetra_hubs) );
		CUDA_SAFE_CALL( cudaFree(_d_tetra_hub_neighs) );
	}
}

template<typename number, typename number4>
void CUDAStarrInteraction<number, number4>::get_settings(input_file &inp) {
	StarrInteraction<number>::get_settings(inp);

	int sort_every;
	if(getInputInt(&inp, "CUDA_sort_every", &sort_every, 0) == KEY_FOUND) {
		if(sort_every > 0) throw oxDNAException("Starr interaction is not compatible with particle sorting, aborting");
	}
}

template<typename number, typename number4>
void CUDAStarrInteraction<number, number4>::cuda_init(number box_side, int N) {
	CUDABaseInteraction<number, number4>::cuda_init(box_side, N);
	StarrInteraction<number>::init();

	_setup_tetra_hubs();

	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N_hubs, &_N_hubs, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N, &N, sizeof(int)) );
	float f_copy = box_side;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_box_side, &f_copy, sizeof(float)) );
	f_copy = this->_fene_sqr_r0;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_fene_sqr_r0, &f_copy, sizeof(float)) );
	f_copy = this->_lin_k;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_lin_k, &f_copy, sizeof(float)) );
	f_copy = this->_fene_K;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_fene_K, &f_copy, sizeof(float)) );
	f_copy = this->_sqr_rcut;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_sqr_rcut, &f_copy, sizeof(float)) );

	COPY_ARRAY_TO_CONSTANT(MD_LJ_sigma, this->_LJ_sigma, 3);
	COPY_ARRAY_TO_CONSTANT(MD_LJ_sqr_sigma, this->_LJ_sqr_sigma, 3);
	COPY_ARRAY_TO_CONSTANT(MD_LJ_sqr_rcut, this->_LJ_sqr_rcut, 3);
	COPY_ARRAY_TO_CONSTANT(MD_LJ_E_cut, this->_LJ_E_cut, 3);
	COPY_ARRAY_TO_CONSTANT(MD_der_LJ_E_cut, this->_der_LJ_E_cut, 3);
}

template<typename number, typename number4>
void CUDAStarrInteraction<number, number4>::_setup_tetra_hubs() {
	BaseParticle<number> **particles = new BaseParticle<number> *[this->_N];
	StarrInteraction<number>::allocate_particles(particles, this->_N);
	int N_strands;
	StarrInteraction<number>::read_topology(this->_N, &N_strands, particles);

	_N_hubs = this->_N_tetramers;
	_h_tetra_hubs = new int[_N_hubs];
	_h_tetra_hub_neighs = new tetra_hub_bonds[_N_hubs];

	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&_d_tetra_hubs, _N_hubs*sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<tetra_hub_bonds>(&_d_tetra_hub_neighs, _N_hubs*sizeof(tetra_hub_bonds)) );

	for(int i = 0; i < this->_N_tetramers; i++) {
		int idx_tetra = this->_N_per_tetramer*i;
		for(int j = 0; j < HUB_SIZE; j++) {
			int idx_hub = idx_tetra + j*this->_N_per_strand;
			CustomParticle<number> *p = static_cast<CustomParticle<number> *>(particles[idx_hub]);
			int rel_idx_hub = i*HUB_SIZE + j;
			_h_tetra_hubs[rel_idx_hub] = p->index;

			// now load all the tetra_hub_bonds structures by looping over all the bonded neighbours
			int nn = 0;
			for(typename set<CustomParticle<number> *>::iterator it = p->bonded_neighs.begin(); it != p->bonded_neighs.end(); it++, nn++) {
				_h_tetra_hub_neighs[rel_idx_hub].n[nn] = (*it)->index;
			}
		}
	}

	CUDA_SAFE_CALL( cudaMemcpy(_d_tetra_hubs, _h_tetra_hubs, _N_hubs*sizeof(int), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_tetra_hub_neighs, _h_tetra_hub_neighs, _N_hubs*sizeof(tetra_hub_bonds), cudaMemcpyHostToDevice) );

	for(int i = 0; i < this->_N; i++) delete particles[i];
	delete[] particles;
}

template<typename number, typename number4>
void CUDAStarrInteraction<number, number4>::compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds) {
	CUDASimpleVerletList<number, number4> *_v_lists = dynamic_cast<CUDASimpleVerletList<number, number4> *>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) throw oxDNAException("use_edge unsupported by StarrInteraction");
		else {
			Starr_forces<number, number4>
				<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
				(d_poss, d_forces, _v_lists->_d_matrix_neighs, _v_lists->_d_number_neighs, d_bonds);
			CUT_CHECK_ERROR("forces_second_step simple_lists error");
		}
	}
	else {
		CUDANoList<number, number4> *_no_lists = dynamic_cast<CUDANoList<number, number4> *>(lists);

		if(_no_lists != NULL) {
			Starr_forces<number, number4>
				<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
				(d_poss,  d_forces, d_bonds);
			CUT_CHECK_ERROR("forces_second_step no_lists error");
		}
	}

	tetra_hub_forces<number, number4>
		<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
		(d_poss, d_forces, _d_tetra_hubs, _d_tetra_hub_neighs);
	CUT_CHECK_ERROR("forces_second_step simple_lists error");
}

extern "C" IBaseInteraction<float> *make_CUDAStarrInteraction_float() {
	return new CUDAStarrInteraction<float, float4>();
}

extern "C" IBaseInteraction<double> *make_CUDAStarrInteraction_double() {
	return new CUDAStarrInteraction<double, LR_double4>();
}

template class CUDAStarrInteraction<float, float4>;
template class CUDAStarrInteraction<double, LR_double4>;