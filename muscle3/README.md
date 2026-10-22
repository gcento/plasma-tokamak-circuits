# MUSCLE3 actor: linear plant models

`muscle_linear_model_actor.m` runs one model of `../linear_models/` as a MUSCLE3 actor
(libmuscle 0.10, Python in-process from MATLAB) with the plant code of `../plant/`, in
place of NICE direct evolutive (`nice_evo_rd`) in controller tests. In IMAS-PDS it is the
program `linear_model` (`workflows/lib/local_programs.ymmsl`) of the workflow
`linear_plant_controller`, run from `local_installs/plasmaless-tokamak-circuits`
(installed by `setup_files/setup_plasmaless.sh`).

## Ports

| Operator | Port | Content |
|---|---|---|
| F_INIT | `equilibrium_f_init` | multi-time equilibrium (reference trajectory); its message timestamp is the start time t0 |
| F_INIT | `pf_active_f_init` | multi-time pf_active, `homogeneous_time` 1, coil currents and voltages |
| S | `pf_active_s` | one-slice pf_active with the 14 power-supply voltages at t, `coil{k}.voltage` = supply k in the order `CS3U CS2U CS1U CS1L CS2L CS3L PF1 .. PF5 VS1 PF6 VS3` (coil names not needed) |
| O_I | `equilibrium_o_i` | one-slice equilibrium at t (below) |
| O_I | `pf_active_o_i` | one-slice coil currents at t (and the voltages applied over the step) |
| O_I | `pf_passive_o_i` | one-slice passive currents at t, one loop per passive state of the model |

Flag `KEEPS_NO_STATE_FOR_NEXT_USE`.

## Settings

| Setting | | Content |
|---|---|---|
| `dt` | mandatory | model time step [s] (zero-order hold) |
| `t_interval` | mandatory | exchange period [s], must equal `dt` |
| `t_end` | mandatory | last output time [s], t0 < `t_end` <= last time of the F_INIT equilibrium |
| `linear_model_dir` | mandatory | absolute path of a model folder of `linear_models/` |
| `md_uri` | optional | IMAS URI of the DD4 machine description (coil names and order); default `data/md_dd4` |
| `operating_point` | optional | `descriptor` (default) or `f_init`: which source wins for the operating point |

A wrong value stops the actor with `linear_model_actor:settings`.

## Operating point and initial state

A deviation model starts at t0 with x = 0 and takes its operating point from two sources:
the descriptor (`operating_point` of `linear_model.json`) and the F_INIT data at t0 (coil
currents of the pf_active sample closest to t0, outputs of the equilibrium slice closest
to t0).
- `descriptor`: the descriptor wins, F_INIT fills the rest. The coil currents and Ip start
  at the linearisation point, so the controller sees a step error at t0.
- `f_init`: F_INIT wins, and the outputs start on the scenario.

The power-supply voltages come from the descriptor only (zero when absent). Each step uses
du = V − op.voltage, so the operating point is held until the controller changes the
voltages.

An absolute model (`plasmaless_model`) has no operating point. It starts from the F_INIT
coil currents at t0 projected on the circuit states (`plant/linear_initial_state.m`), with
passive currents 0, and ignores `operating_point`. Currents with CS1U ≠ CS1L or
VS3U ≠ −VS3L (beyond 1e-6 of the largest current) stop it with `linear_model:circuits`.

Outputs:
- `pf_active_o_i`: op.coil_current + deviations;
- `pf_passive_o_i`: the passive states, named as in the descriptor;
- `equilibrium_o_i`: the F_INIT `ip` and `boundary.geometric_axis` interpolated at t
  (reference pass-through), plus the outputs known to the model
  (`plant/plant_observables_to_equilibrium.m`): Ip replaces `ip`; magnetic axis,
  X-point and gaps are added. An output known by neither source is not sent.

Log of the instance (`stdout.txt`): one line for the model (folder, states, t0, t_end,
dt), one for the operating point (source, coil currents from F_INIT, outputs sent and not
sent), one at the end (number of exchanges). One more line appears only if the F_INIT
pf_active has no sample at t0.

## Timing

As `nice_imas_evo_rd_muscle3` (NICE 3.0.0.dev446):
- t0 = timestamp of `equilibrium_f_init`.
- The first step applies the operating-point voltages (du = 0).
- Each step sends the outputs at t = t_prev + dt on the three O_I ports, with
  next_timestamp t + dt, or None when t + dt > t_end + 1e-9; then the actor stops.
- Otherwise the actor receives `pf_active_s`, which must answer the measurement just sent
  (timestamp t within 1e-6 s, else `linear_model_actor:time`), and applies its voltages
  over the next step.

## Assumptions

- The F_INIT pf_active and the machine description list the same coils in the same order,
  because the controller works by index (`linear_model_actor:coilOrder` otherwise).
- The passive currents are zero at t0: the scenario data has none.
- After the actor returns, the launcher ends MATLAB with `os._exit(0)`: MATLAB R2025b with
  IMAS-MATLAB 5.6.0 segfaults on exit after reading an `imas:hdf5` entry.

## Standalone test

`actor_standalone.ymmsl` couples the actor to `voltage_driver.py` only (no PDS, no
controller). The driver applies a waveform on one power supply, the others at zero, with
the sampling of `utils/supply_sequence_response.m`. F_INIT is all zero and the actor runs
with `operating_point: f_init`, so every output is a deviation from the operating point of
the model, as in MATLAB.

    srun --partition=all_debug --time=00:20:00 --cpus-per-task=2 --mem=8G \
        muscle3/run_actor_standalone.sh [model_dir] [key=value ...]

- `model_dir`: a model folder (default `linear_models/plasmaless_model`).
- Keys: `channel`, `waveform` (`step`, `pulse`, `bipolar`), `amplitude` [V], `t_on` [s]
  (>= dt), `width` [s], `t_end` [s], `dt` [s]. Defaults in `actor_standalone.ymmsl`:
  PF2, step, 10 V, t_on 0.01 s, width 0.1 s, t_end 0.5 s, dt 1 ms.
- Results: `/scratch/users/$USER/linear_model/m3_<date>/channel_response.mat`, next to the
  MUSCLE3 logs.

To check the actor, run `utils/run_supply_sequence.m` with the same model, that supply and
the same waveform, and set `muscle3_results` to that file. The MUSCLE3 curves (dashed) must
lie on the MATLAB ones, and the console line gives their largest relative difference.
