# Linear plant models for IMAS-PDS

Linear models of the ITER coils, vessel and plasma, and the MUSCLE3 actor that runs them
as the plant of IMAS-PDS controller tests (workflow `linear_plant_controller`), in place
of NICE evolutive. Every model has the interface of NICE: the 14 power-supply voltages in,
the ITER circuits inside, the 14 coil currents, the passive currents and the plasma
outputs out. The controller therefore needs no change when the full NICE replaces a
linear model.

## Folders

```
muscle3/         MUSCLE3 actor (muscle_linear_model_actor.m) and its standalone test
plant/           what the actor runs: model loader (configure_linear_model.m), step
                 (plant_step_ids.m), outputs <-> equilibrium IDS, .npy reader
linear_models/   the models, one folder each (linear_model.json + A, B, C), and their contract
generators/      how the models are made
utils/           launchers (MATLAB scripts with the parameters at the top) and their tools
tests/           unit tests
data/md_dd4/     ITER machine description, DD 4.1.1 (HDF5)
```

IMAS-PDS uses three paths, which must stay stable: `muscle3/muscle_linear_model_actor.m`,
`linear_models/<model>/` and `data/md_dd4/`.

## Quick start (MATLAB on a compute node, repository root)

```matlab
add_linear_model_paths                     % plant/ and utils/
load_machine_description();                    % R2024a: before the first runtests
results = runtests('tests'); table(results)
run utils/run_supply_sequence.m          % after editing its parameters
```

## Open-loop test of the models

`utils/run_supply_sequence.m`, parameters at the top: a voltage on each chosen power
supply in turn (`bipolar`, `pulse` or `step`, `spacing` apart) and the response of the 12
circuit currents, one colour per model.
- Computed as the plant computes it (`utils/supply_sequence_response.m`): same loader,
  same discretisation, from x = 0. Every curve is a deviation from the operating point of
  its model.
- The circuit currents are those the controller measures: CS1 = (CS1U + CS1L)/2,
  VS3 = (VS3U − VS3L)/2 (`utils/circuit_currents.m`).
- Output: one line per model (states, growth rate), one figure with the supply voltage and
  the 12 circuit currents, one figure with zmag for the models with a plasma, and
  optionally a `.mat` file.
- An open-loop unstable model (niceLin, vertical_circuit) stops when a coil current
  exceeds `i_limit`: keep the run short (`t_end`) or use the stable models.

## Test through the MUSCLE3 actor

`muscle3/run_actor_standalone.sh [model] [key=value ...]`, in Slurm, applies a waveform
on one power supply through MUSCLE3 (`muscle3/README.md`). Give its results file as
`muscle3_results` to `run_supply_sequence.m`, with the same model, supply and waveform:
the two results must coincide.

## Building models

`generators/README.md`. Launchers: `utils/make_plasmaless_model.m` (machine description),
`utils/make_vertical_circuit_model.m` (equilibrium, e.g. the `out_nice` of an
`inverse_convergence` run), `utils/make_nice_models.m` (NICE deliveries). Each prints what
the plant sees in the new folder (`utils/model_summary.m`) and, with `compare_with`
(`compare_root` for NICE), checks that it reproduces a stored folder
(`utils/same_model.m`).

## Unit tests

- `tests/test_linear_models.m`: loader and step of the plant, descriptor errors, the
  contract on every folder of `linear_models/` (new folders included), the open-loop tool.
- `tests/test_generators.m`: the generators reproduce the stored plasma-free and NICE
  models.

## PDS integration

- The IMAS-PDS program `linear_model` (`workflows/lib/local_programs.ymmsl`) runs
  `muscle3/muscle_linear_model_actor.m` from `$PDS_REPO/local_installs/plasmaless-tokamak-circuits`
  (`setup_files/setup_plasmaless.sh`, branch `IMAS_Muscle3` of
  github.com/gcento/plasma-tokamak-circuits), or from `$LINEAR_MODEL_REPO` when set.
- `linear_plant_controller` selects the model with `plant.linear_model_dir` in
  `cases/overrides/linear_plant_controller_<shot>.ymmsl` (`linear_models/README.md`).
- The actor reads this checkout at run time: edits here reach existing cases.

## Machine description

`data/md_dd4` (DD 4.1.1, HDF5): `em_coupling` (mutual inductances and coil-to-grid flux
maps on the DINA grid), `pf_active` (14 coils with resistances), `pf_passive` (102 passive
loops) and `wall`. `generators/md/convert_md_dd3_to_dd4.py` made it once from the DD3
sources: DINA `em_coupling`/`pf_active`/`equilibrium` of
`/work/projects/dina/SRO_JINTRAC/15MA_10perc/output`, `pf_passive` and `wall` of ITER_MD
115005/3 and 116000/2 (needs `module load IMAS-MUSCLE3/1.0.0-intel-2025b-pds`).
