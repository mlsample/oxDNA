/**
 * @file    ColvarsManager.h
 * @brief   oxDNA-side facade over the Colvars library (engine-proxy interface).
 *
 * The MD backends (CPU and CUDA) talk to Colvars only through this class; the Colvars-side
 * implementation lives in colvarproxy_oxdna.{h,cpp}. User documentation: docs/source/colvars.md.
 *
 * This header must stay includable from nvcc-compiled translation units (.cu): it includes
 * only oxDNA headers and forward-declares every Colvars type. No Colvars header may be
 * included here.
 *
 * Code that uses this class is wrapped in `#ifdef OXDNA_COLVARS` (defined by CMake when
 * -DCOLVARS=ON). The header itself is not guarded so that it can always be parsed.
 */

#ifndef OXDNA_COLVARS_MANAGER_H_
#define OXDNA_COLVARS_MANAGER_H_

// standard headers first: upstream defs.h is not self-contained (needs <string>, <cstdint>)
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "../defs.h"            // number, llint, LR_vector, input_file, oxDNAException
#include "../Boxes/BaseBox.h"

class BaseParticle;
class colvarproxy_oxdna;   // defined in src/Colvars/colvarproxy_oxdna.h (Colvars C++17 code)

/**
 * @brief Unit conversion constants for the Colvars "oxdna" unit system.
 *
 * Colvars internal units == oxDNA simulation units (SU) for length, energy and force,
 * so positions and forces cross the interface with NO conversion. These constants are what
 * the proxy hands to Colvars so that its unit-aware features (kcal/mol / Angstrom literals,
 * kB*T, fs timestep, Langevin extended-system integrators) are consistent.
 */
namespace oxdna_colvars_units {
	/// 1 SU length = 8.518 Angstrom (0.8518 nm)  -> colvarproxy::angstrom_value_ = 1/8.518
	constexpr double ANGSTROM_PER_SU = 8.518;
	/// 1 SU temperature (T_SU = 0.1 <-> 300 K) = 3000 K, as in Utils::get_temperature()
	constexpr double KELVIN_PER_SU = 3000.0;
	/// molar gas constant in kcal/(mol K) (CODATA 2018, exact: R = N_A k_B)
	constexpr double R_KCAL_MOL_K = 1.987204259e-3;
	/// 1 SU energy = k_B * 3000 K = R * 3000 K = 5.961613 kcal/mol (= 4.1419e-20 J)
	/// -> colvarproxy::kcal_mol_value_ = 1/5.961613 = 0.16773984
	constexpr double KCAL_MOL_PER_SU = R_KCAL_MOL_K * KELVIN_PER_SU;
	/// colvarproxy::boltzmann_ in SU energy / K  (exactly 1/3000; consistent with the above)
	constexpr double BOLTZMANN_SU_PER_K = 1.0 / KELVIN_PER_SU;
	/// 1 SU time = 3.03 ps -> Colvars timestep (fs) = dt_SU * 3030
	constexpr double FS_PER_SU_TIME = 3030.0;
}

/**
 * @brief Owns one Colvars module instance (through colvarproxy_oxdna) for an MD run.
 *
 * Ownership: one instance per MD backend, held as `std::unique_ptr<ColvarsManager> _colvars`
 * by MDBackend (protected), so MD_CPUBackend and MD_CUDABackend share it.
 *
 * Index conventions:
 *   - Colvars atom number = oxDNA particle index + 1 (1-based, as in Colvars config files).
 *   - requested_indices() returns 0-based oxDNA particle indices (BaseParticle::index), in
 *     Colvars atom-slot order. All position/force vectors exchanged with compute() use this
 *     order and have exactly requested_indices().size() entries.
 *
 * Units: everything that crosses this interface is in oxDNA simulation units
 * (positions SU length, forces SU energy / SU length, energies SU, temperature T_SU, dt SU time).
 *
 * Step convention: `step` is the MD step number of the configuration whose positions are
 * passed (see DESIGN.md "Step counter"). The engine's step is authoritative over Colvars'.
 *
 * Errors: every failure (bad config, Colvars error, protocol violation) throws oxDNAException.
 */
class ColvarsManager {
public:
	ColvarsManager();
	~ColvarsManager();   // out-of-line (unique_ptr to incomplete colvarproxy_oxdna)

	ColvarsManager(const ColvarsManager &) = delete;
	ColvarsManager &operator=(const ColvarsManager &) = delete;

	/**
	 * @brief Read the colvars_* input keys and, if `colvars_file` is set, create and configure
	 * the Colvars module.
	 *
	 * Must be called once, from the backend's init(), after the configuration has been read
	 * (CONFIG_INFO->curr_step holds the start step, the box is initialised) and before the first
	 * force computation. If `colvars_file` is absent this is a no-op and enabled() == false.
	 *
	 * Sequence (implemented in ColvarsManager.cpp / colvarproxy_oxdna.cpp):
	 *   parse keys -> proxy->set_default_restart_frequency(colvars_restart_frequency) ->
	 *   new colvarmodule -> set_unit_system("oxdna") -> boltzmann_, set_target_temperature(K),
	 *   set_integration_timestep(fs) -> seed RNG -> set prefixes -> read_config_file(colvars_file)
	 *   -> set_initial_step(CONFIG_INFO->curr_step) -> setup_input() (colvars_input_state) ->
	 *   step reconciliation -> setup_output() -> cache requested_indices().
	 *
	 * @param inp       the simulation input (CONFIG_INFO->sim_input)
	 * @param particles all particles (index i <-> particles[i]->index == i)
	 * @param box       simulation box (used for the initial boundaries)
	 * @param dt        MD timestep in SU time
	 * @param T_su      target temperature in SU (K = T_su * 3000)
	 */
	void init(input_file &inp, std::vector<BaseParticle *> &particles, BaseBox *box, number dt, number T_su);

	/// true iff init() found `colvars_file` and the Colvars module was created
	bool enabled() const;

	/// 0-based particle indices Colvars needs, in Colvars atom-slot order. Valid after init().
	const std::vector<int> &requested_indices() const;

	/**
	 * @brief Run one Colvars evaluation and return the bias forces.
	 *
	 * @param step          MD step of the configuration in abs_positions
	 * @param abs_positions UNWRAPPED (absolute, continuous across fix_diffusion) positions of
	 *                      requested_indices(), same order, SU length
	 * @param box           current box; its box_sides() become Colvars' orthorhombic
	 *                      boundaries for this step (Colvars computes minimum images itself)
	 * @param forces_out    resized to requested_indices().size(); receives the Colvars applied
	 *                      forces (SU), same order. The caller ADDS them to the particle forces.
	 * @return bias energy of this step (SU)
	 *
	 * Step rules:
	 *   - first call after init(): evaluates at the initial step;
	 *   - step == previous step + 1: advances Colvars' counter (cvmodule->it++) and evaluates;
	 *   - step == previous step: does NOT call Colvars again; returns the cached forces/energy
	 *     of that step (idempotent re-evaluation, e.g. a backend that recomputes forces);
	 *   - anything else: throws oxDNAException.
	 */
	number compute(llint step, const std::vector<LR_vector> &abs_positions, const BaseBox *box, std::vector<LR_vector> &forces_out);

	/**
	 * @brief Update the Colvars target temperature (e.g. on the CONFIG_INFO "T_updated" event).
	 * @param T_su temperature in SU; forwarded as T_su * 3000 K.
	 */
	void set_temperature(number T_su);

	/// bias energy (SU) returned by the most recent compute(); 0 before the first call
	number bias_energy() const;

	/// step of the most recent compute(); -1 before the first call
	llint last_step() const;

	/**
	 * @brief Write the Colvars state that goes with the last configuration: proxy->post_run()
	 * writes <colvars_output_prefix>.colvars.state and the bias output files, then flushes all
	 * Colvars streams.
	 *
	 * Called whenever the last configuration is printed and the run may stop there
	 * (SimManager::run, oxpy's print_configuration()), so that the state and lastconf describe
	 * the same step. It can be called repeatedly: it writes only when a step has been evaluated
	 * since the previous write (so a run_steps() after run_complete() is written again).
	 *
	 * @param step step of the printed configuration (should equal last_step(), otherwise a
	 *             warning is logged). Safe to call when !enabled() (no-op).
	 */
	void on_run_end(llint step);

	/**
	 * @brief Default seed of the Colvars RNG when colvars_seed is not set: a fixed hash
	 * (splitmix64) of oxDNA's `seed` input value and of the initial step. It does not touch
	 * oxDNA's drand48 stream, so adding Colvars does not change the oxDNA random numbers.
	 */
	static uint64_t derive_seed(long long oxdna_seed, llint initial_step);

private:
	std::unique_ptr<colvarproxy_oxdna> _proxy;   // owns the colvarmodule
	std::vector<int> _requested_indices;          // 0-based, slot order
	std::vector<LR_vector> _cached_forces;        // forces of _last_step (for repeated-step calls)
	number _last_energy = 0.;
	llint _last_step = -1;
	llint _initial_step = 0;                      // CONFIG_INFO->curr_step at init (S0)
	std::shared_ptr<bool> _alive;                 // guards the CONFIG_INFO "T_updated" callback
	bool _enabled = false;
	llint _finalised_step = -1;                   // step of the last on_run_end() write; -1 = none

	// parsed input keys (see DESIGN.md "Input keys")
	std::string _config_file;           // colvars_file (required to enable)
	std::string _output_prefix;         // colvars_output_prefix, default "colvars"
	std::string _input_state;           // colvars_input_state, default "" (none)
	int _restart_frequency = 0;         // colvars_restart_frequency (> 0), default print_conf_interval
	uint64_t _seed = 0;                 // colvars_seed, default derive_seed(seed, initial step)
};

#endif /* OXDNA_COLVARS_MANAGER_H_ */
