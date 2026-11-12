# Collective variables with Colvars

oxDNA can be built with the [Colvars](https://colvars.github.io) library, which lets MD
simulations use Colvars' collective variables and biasing methods: restraints and walls,
(well-tempered) metadynamics, OPES, extended-Lagrangian ABF (eABF) with the CZAR estimator,
steered MD and more. oxDNA talks to Colvars through an engine "proxy", as NAMD, LAMMPS and GROMACS
do, so the Colvars configuration syntax is the standard one documented in the
[Colvars reference manual](https://colvars.github.io/colvars-refman-lammps/colvars-refman-lammps.html).

## Compiling

```bash
cmake .. -DCOLVARS=ON            # add -DCUDA=ON for the GPU backend
make -j
```

The Colvars sources (and the Lepton expression parser it uses for `customFunction`) are shipped in
`src/extern/colvars`. `-DCOLVARS=ON` requires a C++17 compiler. Without `-DCOLVARS=ON` no Colvars
code is compiled, and an input that sets `colvars_file` is refused.

## Input keys

| Key | Meaning | Default |
|---|---|---|
| `colvars_file` | Colvars configuration file. Setting it enables Colvars. | none |
| `colvars_output_prefix` | Prefix of the Colvars output files (`<prefix>.colvars.traj`, `<prefix>.colvars.state`, bias outputs such as `<prefix>.<bias>.czar.pmf`). | `colvars` |
| `colvars_input_state` | Colvars state file (or its prefix) to restart from. | none |
| `colvars_restart_frequency` | Every how many steps the Colvars state is written; must be > 0. The state is also written together with the last configuration at the end of the run. | `print_conf_interval` |
| `colvars_seed` | Seed of the Colvars random number generator (used e.g. by the extended-Lagrangian thermostat). | derived from `seed` and the initial step; it does not consume oxDNA random numbers |

## Conventions

- **Atoms are nucleotides.** Colvars atom number = oxDNA particle index + 1. The position Colvars
  sees is the nucleotide centre of mass; Colvars forces are applied to it (no torques).
- **Units.** Colvars works in oxDNA simulation units (`units oxdna`, set automatically): lengths,
  energies and forces are passed unscaled. Temperatures are in K (T[K] = T[SU] × 3000, with
  k_B = 1/3000 SU/K), and times in fs (1 SU of time ≈ 3030 fs, so `timestep` = `dt` × 3030). Keep
  this in mind for time-based Colvars keys such as `extendedTimeConstant` and
  `extendedLangevinDamping`.
- **Periodic boundaries.** Colvars applies the minimum-image convention itself, using the oxDNA box.
- **Steps.** Colvars follows oxDNA's step counter, also after a restart (the state file's step is
  overridden, with a warning, if they differ).

## Restarting

Restart from the last configuration and the Colvars state of the previous run:

```
conf_file = last_conf.dat
restart_step_counter = 0
colvars_input_state = colvars        # the previous colvars_output_prefix
```

## Supported and unsupported features

- MD only (CPU and CUDA backends). MC/VMMC/FFS simulations are refused.
- Not supported: the barostat, Lees-Edwards boxes, `equilibration_steps > 0`, and anything that needs
  the total force on the atoms (classic ABF, `outputTotalForce`). Extended-Lagrangian ABF (eABF)
  does not need total forces and is supported.
- CUDA: requires `CUDA_sort_every = 0`. Positions of the Colvars atoms are copied to the host and
  the Colvars forces back to the GPU at every step, which costs a roughly fixed time per step
  (about 0.15 ms on a laptop GPU): significant for small systems, minor for large ones.
- The Colvars bias energy is not included in `energy.dat`; use the bias `outputEnergy` option to
  write it to the Colvars trajectory.

## Notes for eABF and meta-eABF users

- When ABF and metadynamics act on the same extended colvar, set `subtractAppliedForce on` in the
  colvar, so that the ABF estimator samples only the coupling-spring force. Without it the
  metadynamics force enters the ABF force estimate and the free-energy profile is biased.
- To keep the extended coordinate inside the grid (instead of adding a wall on the atoms), use
  `reflectingLowerBoundary on` / `reflectingUpperBoundary on`.
- Free-energy profiles: the CZAR estimate is written to `<prefix>.<abf bias>.czar.pmf`.

## Third-party software

- Colvars is distributed under the GNU Lesser General Public License v3
  (`src/extern/colvars/COPYING.LESSER`). Please cite it when you use it
  (see the Colvars documentation for the current reference).
- Lepton is distributed under the MIT license (`src/extern/colvars/lepton/LICENSE`).
- The exact versions shipped are listed in `src/extern/colvars/VERSION`.
