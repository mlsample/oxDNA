/**
 * @file    colvars_kernels.cuh
 * @brief   Device kernels for the Colvars host round trip of the CUDA MD backends.
 *
 * Included by MD_CUDABackend.cu and MD_CUDAMixedBackend.cu. All kernels are templates, so the
 * header can be included by more than one translation unit without duplicate symbols.
 *
 * Index convention: idx[k] is the 0-based oxDNA particle index of Colvars slot k. The CUDA
 * backends refuse Colvars with CUDA_sort_every > 0, so the GPU slot of a particle equals its
 * index and idx[k] addresses the device arrays directly.
 */

#ifndef COLVARS_KERNELS_CUH_
#define COLVARS_KERNELS_CUH_

#include "../cuda_defs.h"

/// threads per block for the (tiny) Colvars gather / scatter launches
#define COLVARS_KERNEL_BLOCK 128

/**
 * @brief Copy the xyz components of the requested particles into a slot-ordered double buffer.
 *
 * T is c_number4 (float4 or LR_double4) for MD_CUDABackend, LR_double4 for the mixed backend.
 * The source is only read; the .w component (packed btype/index) is never touched.
 */
template<typename T>
__global__ void colvars_gather_positions(const T *poss, const int *idx, int n, LR_double4 *out) {
	const int k = blockIdx.x * blockDim.x + threadIdx.x;
	if(k >= n) return;

	const T p = poss[idx[k]];
	LR_double4 r;
	r.x = (double) p.x;
	r.y = (double) p.y;
	r.z = (double) p.z;
	r.w = 0.f;
	out[k] = r;
}

/**
 * @brief Add the slot-ordered Colvars forces to the particle force array (xyz only).
 *
 * .w of forces carries the per-particle energy used by CUDA_print_energy and is left untouched.
 * Colvars slots refer to distinct particles, so there are no write conflicts.
 */
template<typename F4>
__global__ void colvars_add_forces(F4 *forces, const int *idx, int n, const F4 *colvars_forces) {
	const int k = blockIdx.x * blockDim.x + threadIdx.x;
	if(k >= n) return;

	const int i = idx[k];
	const F4 f = colvars_forces[k];
	F4 F = forces[i];
	F.x += f.x;
	F.y += f.y;
	F.z += f.z;
	forces[i] = F;
}

#endif /* COLVARS_KERNELS_CUH_ */
