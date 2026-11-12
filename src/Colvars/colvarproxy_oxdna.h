/**
 * @file    colvarproxy_oxdna.h
 * @brief   Colvars engine proxy for oxDNA (implementation of colvarproxy).
 *
 * Built only with -DCOLVARS=ON. This header includes Colvars headers (C++17) and must NOT be
 * included from nvcc-compiled translation units: the backends talk to ColvarsManager only.
 *
 * Conventions (see docs/colvars/DESIGN.md):
 *   - Colvars atom number = oxDNA particle index + 1; the internal atom id is the 0-based index.
 *   - Colvars internal units == oxDNA simulation units ("oxdna" unit system).
 *   - Positions handed in are unwrapped; Colvars applies its own minimum image through
 *     boundaries_ (use_internal_pbc_ = true), refreshed from the oxDNA box at every evaluation.
 *   - Every Colvars error throws oxDNAException.
 */

#ifndef OXDNA_COLVARPROXY_OXDNA_H_
#define OXDNA_COLVARPROXY_OXDNA_H_

#include <cstdint>
#include <random>
#include <string>
#include <vector>

// oxDNA first: Colvars (a system include) then silently redefines PI for the rest of this TU
#include "../defs.h"

#include "colvarmodule.h"
#include "colvarproxy.h"

class colvarproxy_oxdna: public colvarproxy {
public:
	/// Everything the proxy needs to create and configure the module.
	struct settings {
		std::string config_file;        ///< Colvars configuration file
		std::string output_prefix;      ///< output and restart prefix
		std::string input_state;        ///< prefix (or path) of the state file to read; may be empty
		int restart_frequency = 0;      ///< default restart frequency (steps)
		uint64_t seed = 0;              ///< seed of the Gaussian RNG
		double temperature_K = 0.;      ///< target temperature (K)
		double timestep_fs = 1.;        ///< integration timestep (fs)
		long long initial_step = 0;     ///< oxDNA step of the first evaluation
		int n_particles = 0;            ///< number of particles (valid atom numbers are 1..N)
	};

	colvarproxy_oxdna();
	~colvarproxy_oxdna() override;

	colvarproxy_oxdna(const colvarproxy_oxdna &) = delete;
	colvarproxy_oxdna &operator=(const colvarproxy_oxdna &) = delete;

	/// Create the colvarmodule, read the configuration and the state file, open the outputs.
	void setup_module(const settings &s);

	/**
	 * @brief One Colvars evaluation.
	 *
	 * @param abs_positions unwrapped positions, in atom-slot order (SU)
	 * @param box_sides     orthorhombic box sides (SU)
	 * @param advance       if true, increment the Colvars step counter before evaluating
	 * @param forces_out    applied forces, in atom-slot order (SU)
	 * @return the bias energy (SU)
	 */
	double compute(const std::vector<LR_vector> &abs_positions, const LR_vector &box_sides, bool advance, std::vector<LR_vector> &forces_out);

	/// Write the final state and the bias output files, flush all streams.
	void finalise();

	/// 0-based particle indices, in slot order
	const std::vector<int> &atom_ids() const {
		return atoms_ids;
	}

	cvm::step_number current_step() const;

	void update_temperature_K(double T_K);

	// colvarproxy interface
	int set_unit_system(std::string const &units_in, bool check_only) override;
	void request_total_force(bool yesno) override;
	bool total_forces_enabled() const override {
		return false;
	}
	bool total_forces_same_step() const override {
		return false;
	}
	void log(std::string const &message) override;
	void error(std::string const &message) override;
	cvm::real rand_gaussian() override;
	void add_energy(cvm::real energy) override {
		_bias_energy += energy;
	}

	int init_atom(int atom_number) override;
	int init_atom(cvm::residue_id const &residue, std::string const &atom_name, std::string const &segment_id) override {
		return colvarproxy::init_atom(residue, atom_name, segment_id);
	}
	int check_atom_id(int atom_number) override;
	int check_atom_id(cvm::residue_id const &residue, std::string const &atom_name, std::string const &segment_id) override {
		return colvarproxy::check_atom_id(residue, atom_name, segment_id);
	}

private:
	int _n_particles = 0;
	double _bias_energy = 0.;
	/// errors throw only while this is true (never from a destructor)
	bool _throw_on_error = true;
	std::mt19937_64 _rng;
	std::normal_distribution<double> _gauss {0., 1.};
};

#endif /* OXDNA_COLVARPROXY_OXDNA_H_ */
