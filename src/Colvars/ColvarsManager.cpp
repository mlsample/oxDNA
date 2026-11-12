/**
 * @file    ColvarsManager.cpp
 * @brief   Implementation of the oxDNA-side facade over Colvars (see ColvarsManager.h).
 */

#include "ColvarsManager.h"

#include <cstdlib>

#include "colvarproxy_oxdna.h"

#include "../Particles/BaseParticle.h"
#include "../Utilities/ConfigInfo.h"
#include "../Utilities/Logger.h"
#include "../Utilities/oxDNAException.h"
#include "../Utilities/parse_input/parse_input.h"

ColvarsManager::ColvarsManager() {

}

ColvarsManager::~ColvarsManager() {
	if(_alive) {
		*_alive = false;
	}
}

void ColvarsManager::init(input_file &inp, std::vector<BaseParticle *> &particles, BaseBox *box, number dt, number T_su) {
	if(_enabled) {
		throw oxDNAException("Colvars: ColvarsManager::init() called twice");
	}

	if(getInputString(&inp, "colvars_file", _config_file, 0) == KEY_NOT_FOUND || _config_file.empty()) {
		return;
	}

	_output_prefix = "colvars";
	getInputString(&inp, "colvars_output_prefix", _output_prefix, 0);
	if(_output_prefix.empty()) {
		throw oxDNAException("Colvars: colvars_output_prefix cannot be empty");
	}
	getInputString(&inp, "colvars_input_state", _input_state, 0);

	long long restart_freq = 0;
	if(getInputLLInt(&inp, "colvars_restart_frequency", &restart_freq, 0) == KEY_NOT_FOUND) {
		getInputLLInt(&inp, "print_conf_interval", &restart_freq, 0);
	}
	if(restart_freq < 0) {
		throw oxDNAException("Colvars: colvars_restart_frequency must be >= 0 (0 = write the state only at the end)");
	}
	_restart_frequency = (int) restart_freq;

	if(getInputLLInt(&inp, "colvars_seed", &_seed, 0) == KEY_NOT_FOUND) {
		_seed = (long long) lrand48();
	}

	_initial_step = CONFIG_INFO->curr_step;

	colvarproxy_oxdna::settings s;
	s.config_file = _config_file;
	s.output_prefix = _output_prefix;
	s.input_state = _input_state;
	s.restart_frequency = _restart_frequency;
	s.seed = (uint64_t) _seed;
	s.temperature_K = (double) T_su * oxdna_colvars_units::KELVIN_PER_SU;
	s.timestep_fs = (double) dt * oxdna_colvars_units::FS_PER_SU_TIME;
	s.initial_step = _initial_step;
	s.n_particles = (int) particles.size();

	for(size_t i = 0; i < particles.size(); i++) {
		if(particles[i]->index != (int) i) {
			throw oxDNAException("Colvars: particle %zu has index %d; Colvars needs particles[i]->index == i", i, particles[i]->index);
		}
	}

	OX_LOG(Logger::LOG_INFO, "Colvars: config = '%s', output prefix = '%s', input state = '%s', restart frequency = %d, seed = %lld, T = %g K, dt = %g fs, initial step = %lld",
			_config_file.c_str(), _output_prefix.c_str(), _input_state.c_str(), _restart_frequency, _seed, s.temperature_K, s.timestep_fs, _initial_step);

	_proxy = std::make_unique<colvarproxy_oxdna>();
	_proxy->setup_module(s);

	_requested_indices = _proxy->atom_ids();
	for(int idx : _requested_indices) {
		if(idx < 0 || idx >= (int) particles.size()) {
			throw oxDNAException("Colvars: requested particle index %d is out of range", idx);
		}
	}
	_cached_forces.assign(_requested_indices.size(), LR_vector(0., 0., 0.));
	_last_energy = 0.;
	_last_step = -1;
	_finalised = false;
	_enabled = true;

	// follow temperature changes (e.g. T updated from oxpy)
	_alive = std::make_shared<bool>(true);
	std::shared_ptr<bool> alive = _alive;
	CONFIG_INFO->subscribe("T_updated", [this, alive]() {
		if(*alive && _enabled) {
			this->set_temperature(CONFIG_INFO->temperature());
		}
	});

	OX_LOG(Logger::LOG_INFO, "Colvars: initialised, %zu particles requested", _requested_indices.size());
}

bool ColvarsManager::enabled() const {
	return _enabled;
}

const std::vector<int> &ColvarsManager::requested_indices() const {
	return _requested_indices;
}

number ColvarsManager::compute(llint step, const std::vector<LR_vector> &abs_positions, const BaseBox *box, std::vector<LR_vector> &forces_out) {
	if(!_enabled) {
		throw oxDNAException("Colvars: compute() called but Colvars is not enabled");
	}

	bool advance;
	if(_last_step < 0) {
		if(step != _initial_step) {
			throw oxDNAException("Colvars: the first evaluation must be at the initial step %lld, not at step %lld", _initial_step, step);
		}
		advance = false;
	}
	else if(step == _last_step) {
		// repeated evaluation of the same configuration: do not sample/deposit twice
		forces_out = _cached_forces;
		return _last_energy;
	}
	else if(step == _last_step + 1) {
		advance = true;
	}
	else {
		throw oxDNAException("Colvars: non-sequential step %lld (the last evaluated step was %lld)", step, _last_step);
	}

	_last_energy = (number) _proxy->compute(abs_positions, box->box_sides(), advance, _cached_forces);
	_last_step = step;

	if((llint) _proxy->current_step() != step) {
		throw oxDNAException("Colvars: internal step counter (%lld) out of sync with oxDNA (%lld)", (long long) _proxy->current_step(), step);
	}

	forces_out = _cached_forces;
	return _last_energy;
}

void ColvarsManager::set_temperature(number T_su) {
	if(!_enabled) {
		return;
	}
	_proxy->update_temperature_K((double) T_su * oxdna_colvars_units::KELVIN_PER_SU);
}

number ColvarsManager::bias_energy() const {
	return _last_energy;
}

llint ColvarsManager::last_step() const {
	return _last_step;
}

void ColvarsManager::on_run_end(llint step) {
	if(!_enabled || _finalised) {
		return;
	}
	if(step != _last_step) {
		OX_LOG(Logger::LOG_WARNING, "Colvars: run ended at step %lld but the last Colvars evaluation was at step %lld", step, _last_step);
	}
	_finalised = true;
	_proxy->finalise();
	OX_LOG(Logger::LOG_INFO, "Colvars: final state written to %s.colvars.state", _output_prefix.c_str());
}
