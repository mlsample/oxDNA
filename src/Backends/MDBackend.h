/**
 * @file    MDBackend.h
 * @date    23/nov/2010
 * @author  lorenzo
 *
 *
 */

#ifndef MDBACKEND_H_
#define MDBACKEND_H_

#include "SimBackend.h"

#ifdef OXDNA_COLVARS
#include <memory>
#include "../Colvars/ColvarsManager.h"
#endif

using namespace std;

/**
 * @brief Abstract class for MD simulations.
 *
 * This class sets up a basic MD simulation. It does not do any sensible computation but
 * takes care of the most basic input/output operations associated with MD simulations.
 *
 * @verbatim
 [reset_initial_com_momentum = <bool> (if true the momentum of the centre of mass of the initial configuration will be set to 0. Defaults to false to enforce the reproducibility of the trajectory)]
 [reset_com_momentum = <bool> (if true the momentum of the centre of mass will be set to 0 each time fix_diffusion is performed. Defaults to false to enforce the reproducibility of the trajectory)]
 [use_barostat = <bool> (apply an MC-like barostat to the simulation)]
 [barostat_probability = <float> (the probability that, at each step, the barostat gets activated)]
 [P = <float> (the pressure of the simulation)]
 [delta_L = <float> (controls the box side change by the MC-like barostat)]
 [colvars_file = <path> (Colvars configuration file; enables Colvars biases on MD runs. Requires oxDNA compiled with -DCOLVARS=ON. Colvars atom numbers are particle index + 1 and Colvars works in oxDNA units ("units oxdna"))]
 [colvars_output_prefix = <string> (prefix of the Colvars output files, e.g. <prefix>.colvars.traj and <prefix>.colvars.state. Defaults to colvars)]
 [colvars_input_state = <path> (Colvars state file (or its prefix) to restart from. Defaults to none)]
 [colvars_restart_frequency = <int> (every how many steps the Colvars state file is written; must be > 0. The state is also written together with the last configuration at the end of the run. Defaults to print_conf_interval)]
 [colvars_seed = <int> (seed of the Colvars Gaussian random number generator. Defaults to a fixed hash of the oxDNA seed and of the initial step, which does not consume oxDNA random numbers)]
 @endverbatim
 */

class MDBackend: public SimBackend {
protected:
	number _dt = 0.;
	bool _refresh_velocities = false;
	bool _reset_initial_com_momentum = false;
	bool _reset_com_momentum = false;

	// Stuff needed by the barostat
	bool _use_barostat = false;
	bool _barostat_isotropic = true;
	bool _barostat_molecular = false;
	number _delta_L = 0.;
	number _barostat_probability = 1.;
	number _barostat_acceptance = 0.;

	// timers
	std::shared_ptr<Timer> _timer_first_step;
	std::shared_ptr<Timer> _timer_forces;
	std::shared_ptr<Timer> _timer_lists;
	std::shared_ptr<Timer> _timer_thermostat;
	std::shared_ptr<Timer> _timer_barostat;

	bool _lees_edwards;
	number _shear_rate;

#ifdef OXDNA_COLVARS
	/// Colvars biases (created in get_settings, initialised by the concrete backend's init())
	std::unique_ptr<ColvarsManager> _colvars;
#endif

	/**
	 * @brief Set the momentum of the centre of mass to 0.
	 */
	virtual void _reset_momentum();

	void _generate_vel();
	virtual bool _is_barostat_active();

public:
	MDBackend();
	virtual ~MDBackend();

	void get_settings(input_file &inp);
	void init();
	void fix_diffusion();

	virtual void print_observables();

#ifdef OXDNA_COLVARS
	void on_run_end() override;
#endif
};

#endif /* MDBACKEND_H_ */
