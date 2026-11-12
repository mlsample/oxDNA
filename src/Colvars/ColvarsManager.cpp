/**
 * @file    ColvarsManager.cpp
 * @brief   Implementation of the oxDNA-side facade over Colvars (see ColvarsManager.h).
 */

#include "ColvarsManager.h"

#include <random>

#include "colvarproxy_oxdna.h"

#include "../Particles/BaseParticle.h"
#include "../Utilities/ConfigInfo.h"
#include "../Utilities/Logger.h"
#include "../Utilities/oxDNAException.h"
#include "../Utilities/parse_input/parse_input.h"
#include "../Utilities/Utils.h"

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
	// Colvars biases (e.g. OPES) compute "step % restart frequency" without checking for zero
	if(restart_freq <= 0) {
		throw oxDNAException("Colvars: colvars_restart_frequency (or print_conf_interval, its default) must be a positive number of steps, not %lld", restart_freq);
	}
	if(restart_freq > INT32_MAX) {
		throw oxDNAException("Colvars: colvars_restart_frequency (or print_conf_interval, its default) must be at most %d, not %lld", INT32_MAX, restart_freq);
	}
	_restart_frequency = (int) restart_freq;

	_initial_step = CONFIG_INFO->curr_step;

	// the default seed must not consume oxDNA's drand48 stream: Colvars on/off would change the trajectory
	long long seed_in = 0;
	std::string seed_source;
	if(getInputLLInt(&inp, "colvars_seed", &seed_in, 0) == KEY_FOUND) {
		_seed = (uint64_t) seed_in;
		seed_source = "colvars_seed";
	}
	else {
		int oxdna_seed = 0;
		if(getInputInt(&inp, "seed", &oxdna_seed, 0) == KEY_FOUND) {
			_seed = derive_seed(oxdna_seed, _initial_step);
			seed_source = Utils::sformat("derived from seed = %d and initial step %lld", oxdna_seed, _initial_step);
		}
		else {
			// oxDNA itself is not reproducible without "seed"
			std::random_device rd;
			_seed = ((uint64_t) rd() << 32) ^ (uint64_t) rd();
			seed_source = "drawn from std::random_device, since seed is not set";
		}
	}

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

	OX_LOG(Logger::LOG_INFO, "Colvars: config = '%s', output prefix = '%s', input state = '%s', restart frequency = %d, T = %g K, dt = %g fs, initial step = %lld",
			_config_file.c_str(), _output_prefix.c_str(), _input_state.c_str(), _restart_frequency, s.temperature_K, s.timestep_fs, _initial_step);
	OX_LOG(Logger::LOG_INFO, "Colvars: random number generator seed = %llu (%s)", (unsigned long long) _seed, seed_source.c_str());

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
	_finalised_step = -1;
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
	if(!_enabled) {
		return;
	}
	if(_last_step < 0) {
		OX_LOG(Logger::LOG_INFO, "Colvars: no step has been evaluated, the Colvars state is not written");
		return;
	}
	if(_last_step == _finalised_step) {
		// nothing new since the last write
		return;
	}
	if(step != _last_step) {
		OX_LOG(Logger::LOG_WARNING, "Colvars: the configuration was printed at step %lld but the last Colvars evaluation was at step %lld", step, _last_step);
	}
	_finalised_step = _last_step;
	_proxy->finalise();
	OX_LOG(Logger::LOG_INFO, "Colvars: state at step %lld written to %s.colvars.state", _last_step, _output_prefix.c_str());
}

uint64_t ColvarsManager::derive_seed(long long oxdna_seed, llint initial_step) {
	// splitmix64 finaliser (Steele, Lea and Flood 2014): a fixed, platform-independent mix
	auto mix = [](uint64_t z) {
		z += 0x9E3779B97F4A7C15ULL;
		z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
		z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
		return z ^ (z >> 31);
	};
	// the constant separates this stream from any other use of the oxDNA seed
	const uint64_t COLVARS_SEED_SALT = 0x436F6C7661727331ULL;   // "Colvars1"
	return mix(mix((uint64_t) oxdna_seed ^ COLVARS_SEED_SALT) ^ (uint64_t) initial_step);
}
