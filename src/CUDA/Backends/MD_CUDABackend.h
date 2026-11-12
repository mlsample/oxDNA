/**
 * @file    MD_CUDABackend.h
 * @date    03/set/2010
 * @author  lorenzo
 *
 *
 */

#ifndef MD_CUDABACKEND_H_
#define MD_CUDABACKEND_H_

#include <cuda.h>
#include <cuda_runtime_api.h>

#include "CUDABaseBackend.h"
#include "../../Backends/MDBackend.h"

#include "../CUDAUtils.h"
#include "../Thermostats/CUDABrownianThermostat.h"
#include "../cuda_utils/cuda_device_utils.h"
#include "../Lists/CUDANoList.h"
#include "../Lists/CUDASimpleVerletList.h"

union CUDA_trap;

/**
 * @brief Manages a MD simulation on GPU with CUDA.
 */
class MD_CUDABackend: public MDBackend, public CUDABaseBackend{
protected:
	bool _use_edge;
	bool _any_rigid_body;
	bool _avoid_cpu_calculations;

	int *_h_gpu_index, *_h_cpu_index;

	std::shared_ptr<Timer> _timer_sorting;

	c_number4 *_d_vels, *_h_vels;
	c_number4 *_d_Ls, *_h_Ls;
	c_number4 *_d_forces, *_h_forces;
	c_number4 *_d_torques, *_h_torques;

	std::vector<int> _h_particles_to_mols;
	int *_d_particles_to_mols, *_d_mol_sizes, *_d_buff_particles_to_mols;
	c_number4 *_d_molecular_coms;

	c_number4 *_d_buff_vels, *_d_buff_Ls;

	llint _barostat_attempts, _barostat_accepted;
	int _update_st_every = 0;
	bool _update_particle_st = false;

	bool _print_energy;

	std::shared_ptr<CUDABaseThermostat> _cuda_thermostat;

	bool _cuda_barostat_always_refresh = false;
	std::shared_ptr<CUDABrownianThermostat> _cuda_barostat_thermostat;

	CUDA_trap *_d_ext_forces;
	int _max_ext_forces;

	virtual void _gpu_to_host();
	virtual void _host_to_gpu();
	virtual void _apply_external_forces_changes();

	virtual void _sort_particles();
	virtual void _rescale_molecular_positions(c_number4 new_Ls, c_number4 old_Ls, bool is_reverse_move);
	virtual void _rescale_positions(c_number4 new_Ls, c_number4 old_Ls);

	virtual void _first_step();
	virtual void _apply_barostat();
	virtual void _forces_second_step();
	virtual void _set_external_forces();

	virtual void _thermalize();
	virtual void _update_stress_tensor();

	virtual void _init_CUDA_MD_symbols();

	/**
	 * Colvars host round trip (used only when oxDNA is built with -DCOLVARS=ON and the input sets
	 * colvars_file; see docs/source/colvars.md). The members are declared
	 * unconditionally so that the class layout does not depend on OXDNA_COLVARS being defined
	 * consistently in every translation unit that includes this header.
	 */
	int _colvars_n = 0;                          ///< number of Colvars atom slots
	int *_d_colvars_idx = nullptr;               ///< device copy of requested_indices()
	LR_double4 *_d_colvars_pos = nullptr;        ///< device gather buffer (slot order)
	c_number4 *_d_colvars_forces = nullptr;      ///< device Colvars forces (slot order)
	LR_double4 *_h_colvars_pos = nullptr;        ///< pinned host mirror of _d_colvars_pos
	c_number4 *_h_colvars_forces = nullptr;      ///< pinned host mirror of _d_colvars_forces
	std::vector<LR_vector> _colvars_offset;      ///< get_abs_pos(p) - p->pos, per slot (unwrapping)
	std::vector<LR_vector> _colvars_abs_pos;     ///< unwrapped positions handed to Colvars
	std::vector<LR_vector> _colvars_f;           ///< forces returned by Colvars
	number _colvars_bias_energy = 0.;            ///< bias energy (SU) of the last evaluation

	/// true iff a Colvars module is active for this run
	bool _colvars_active() const;
	/// creates the Colvars module (if colvars_file is set) and the device buffers
	void _colvars_init();
	/// refreshes _colvars_offset from the host particles (after any host -> GPU upload)
	void _colvars_refresh_offsets();
	/// fills _d_colvars_pos with the positions of the requested particles (overridden by the mixed backend)
	virtual void _colvars_gather_positions();
	/**
	 * Evaluates Colvars for the configuration of MD step `step` and ADDS the bias forces to _d_forces.
	 * Must be called after _set_external_forces() (which overwrites _d_forces). If from_host is true
	 * the positions are taken from the host particles (used once in init()), otherwise they are
	 * gathered from the device (one D2H and one H2D copy).
	 */
	void _colvars_apply(llint step, bool from_host = false);
	void _colvars_free();

public:
	MD_CUDABackend();
	virtual ~MD_CUDABackend();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void sim_step();

	virtual void apply_simulation_data_changes();
	virtual void apply_changes_to_simulation_data();
};

#endif /* MD_CUDABACKEND_H_ */
