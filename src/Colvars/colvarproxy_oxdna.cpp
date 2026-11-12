/**
 * @file    colvarproxy_oxdna.cpp
 * @brief   Colvars engine proxy for oxDNA.
 *
 * Modelled on Colvars' misc_interfaces/stubs/colvarproxy_stub and on colvarproxy_lammps.
 */

#include "colvarproxy_oxdna.h"

#include <sstream>

#include "colvars_version.h"

#include "ColvarsManager.h"
#include "../Utilities/oxDNAException.h"
#include "../Utilities/Logger.h"

colvarproxy_oxdna::colvarproxy_oxdna() {
	engine_name_ = "oxDNA";
	version_int = get_version_from_string(COLVARS_VERSION);

	// oxDNA hands Colvars unwrapped positions and the box: Colvars does the minimum image itself
	use_internal_pbc_ = true;
	b_simulation_running = true;
	b_simulation_continuing = false;
	engine_ready_ = true;

	// "oxdna" unit system: Colvars internal units == oxDNA simulation units
	units = "oxdna";
	angstrom_value_ = 1. / oxdna_colvars_units::ANGSTROM_PER_SU;
	kcal_mol_value_ = 1. / oxdna_colvars_units::KCAL_MOL_PER_SU;
	boltzmann_ = oxdna_colvars_units::BOLTZMANN_SU_PER_K;

	// masses are 1 and charges 0 for every particle, and they never change
	updated_masses_ = updated_charges_ = true;
}

colvarproxy_oxdna::~colvarproxy_oxdna() {
	// the base-class destructor deletes the module: never throw from there
	_throw_on_error = false;
}

void colvarproxy_oxdna::setup_module(const settings &s) {
	if(cvmodule != nullptr) {
		throw oxDNAException("Colvars: the Colvars module has already been initialised");
	}

	_n_particles = s.n_particles;
	_rng.seed(s.seed);
	_gauss.reset();

	// the module constructor copies the default restart frequency: it must be set first
	set_default_restart_frequency(s.restart_frequency);

	cvmodule = new colvarmodule(this);
	cvmodule->log("Using the oxDNA interface (Colvars version " + std::string(COLVARS_VERSION) + ").\n");
	cvmodule->log("oxDNA units: 1 length unit = 8.518 A, 1 energy unit = " + cvm::to_str(oxdna_colvars_units::KCAL_MOL_PER_SU) + " kcal/mol, k_B = 1/3000 per K.\n");

	set_target_temperature(s.temperature_K);
	set_integration_timestep(s.timestep_fs);

	set_output_prefix(s.output_prefix);
	set_restart_output_prefix(s.output_prefix);
	if(!s.input_state.empty()) {
		set_input_prefix(s.input_state);
	}

	if(cvmodule->read_config_file(s.config_file.c_str()) != COLVARS_OK) {
		throw oxDNAException("Colvars: error while reading the configuration file '%s'", s.config_file.c_str());
	}
	if(cvmodule->size() == 0) {
		OX_LOG(Logger::LOG_WARNING, "Colvars: the configuration file '%s' defines no variables or biases", s.config_file.c_str());
	}

	cvmodule->set_initial_step(static_cast<cvm::step_number>(s.initial_step));

	int err = colvarproxy::setup();
	err |= cvmodule->update_engine_parameters();
	err |= cvmodule->setup_input();
	if(err != COLVARS_OK) {
		throw oxDNAException("Colvars: error while setting up the Colvars module (input stage)");
	}

	// oxDNA's step counter is authoritative over the one stored in a state file
	if(cvmodule->it != static_cast<cvm::step_number>(s.initial_step) || cvmodule->it_restart != static_cast<cvm::step_number>(s.initial_step)) {
		OX_LOG(Logger::LOG_WARNING, "Colvars: the state file step (%lld) differs from the oxDNA step (%lld): oxDNA is authoritative", (long long) cvmodule->it, s.initial_step);
		cvmodule->it = cvmodule->it_restart = static_cast<cvm::step_number>(s.initial_step);
	}

	if(cvmodule->setup_output() != COLVARS_OK) {
		throw oxDNAException("Colvars: error while setting up the Colvars module (output stage)");
	}
}

double colvarproxy_oxdna::compute(const std::vector<LR_vector> &abs_positions, const LR_vector &box_sides, bool advance, std::vector<LR_vector> &forces_out) {
	const size_t n = atoms_ids.size();
	if(abs_positions.size() != n) {
		throw oxDNAException("Colvars: %zu positions were passed, but Colvars requested %zu atoms", abs_positions.size(), n);
	}

	if(advance) {
		cvmodule->it++;
	}
	b_simulation_continuing = false;

	boundaries_.set_boundaries(true, true, true, cvm::rvector(box_sides.x, 0., 0.), cvm::rvector(0., box_sides.y, 0.), cvm::rvector(0., 0., box_sides.z));

	for(size_t i = 0; i < n; i++) {
		atoms_positions[i] = cvm::rvector(abs_positions[i].x, abs_positions[i].y, abs_positions[i].z);
		atoms_new_colvar_forces[i].reset();
	}
	_bias_energy = 0.;

	// calc() also runs cvmodule->end_of_step() and colvarproxy::end_of_step()
	if(cvmodule->calc() != COLVARS_OK) {
		throw oxDNAException("Colvars: error in the collective variables module at step %lld", (long long) cvmodule->it);
	}

	forces_out.resize(n);
	for(size_t i = 0; i < n; i++) {
		const cvm::rvector &f = atoms_new_colvar_forces[i];
		forces_out[i] = LR_vector((number) f.x, (number) f.y, (number) f.z);
	}

	return _bias_energy;
}

void colvarproxy_oxdna::finalise() {
	if(cvmodule == nullptr) {
		return;
	}
	if(post_run() != COLVARS_OK) {
		throw oxDNAException("Colvars: error while writing the final Colvars state and output files");
	}
}

cvm::step_number colvarproxy_oxdna::current_step() const {
	return (cvmodule != nullptr) ? cvmodule->it : 0;
}

void colvarproxy_oxdna::update_temperature_K(double T_K) {
	set_target_temperature(T_K);
}

int colvarproxy_oxdna::set_unit_system(std::string const &units_in, bool check_only) {
	// the only unit system is oxDNA's: anything else would silently rescale every number
	if(units_in != "oxdna") {
		error("Error: unit system \"" + units_in + "\" is not supported by the oxDNA interface: only \"units oxdna\" "
				"(oxDNA simulation units: 1 length unit = 0.8518 nm, 1 energy unit = k_B * 3000 K) is allowed.\n");
		return COLVARS_INPUT_ERROR;
	}
	if(!check_only) {
		units = "oxdna";
		angstrom_value_ = 1. / oxdna_colvars_units::ANGSTROM_PER_SU;
		kcal_mol_value_ = 1. / oxdna_colvars_units::KCAL_MOL_PER_SU;
		boltzmann_ = oxdna_colvars_units::BOLTZMANN_SU_PER_K;
	}
	return COLVARS_OK;
}

void colvarproxy_oxdna::request_total_force(bool yesno) {
	if(yesno) {
		error("Error: total forces (e.g. classic ABF, outputTotalForce, TI) are not implemented by the oxDNA interface; "
				"use extended-Lagrangian ABF (eABF), metadynamics or restraints instead.\n");
	}
	total_force_requested = false;
}

void colvarproxy_oxdna::log(std::string const &message) {
	std::istringstream ss(message);
	std::string line;
	while(std::getline(ss, line)) {
		OX_LOG(Logger::LOG_INFO, "Colvars: %s", line.c_str());
	}
}

void colvarproxy_oxdna::error(std::string const &message) {
	add_error_msg(message);
	std::string msg(message);
	while(!msg.empty() && (msg.back() == '\n' || msg.back() == ' ')) {
		msg.pop_back();
	}
	if(_throw_on_error) {
		throw oxDNAException("Colvars: %s", msg.c_str());
	}
	OX_LOG(Logger::LOG_WARNING, "Colvars: %s", msg.c_str());
}

cvm::real colvarproxy_oxdna::rand_gaussian() {
	return _gauss(_rng);
}

int colvarproxy_oxdna::check_atom_id(int atom_number) {
	if(atom_number < 1 || atom_number > _n_particles) {
		error("Error: invalid atom number " + cvm::to_str(atom_number) + ": oxDNA atom numbers are particle index + 1, i.e. 1.." + cvm::to_str(_n_particles) + ".\n");
		return COLVARS_INPUT_ERROR;
	}
	return atom_number - 1;
}

int colvarproxy_oxdna::init_atom(int atom_number) {
	int const aid = check_atom_id(atom_number);
	if(aid < 0) {
		return COLVARS_INPUT_ERROR;
	}

	int index = find_atom_by_id(aid);
	if(index >= 0) {
		atoms_refcount[index] += 1;
		return index;
	}

	index = add_atom_slot(aid);   // mass 1, charge 0
	atoms_masses[index] = 1.;
	atoms_charges[index] = 0.;
	updated_masses_ = updated_charges_ = true;
	return index;
}
