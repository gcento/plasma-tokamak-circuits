# Linear plant models

Linear models of the ITER coils, vessel and plasma. In IMAS-PDS they are selected by the workflow `linear_plant_controller`.

Each folder is one model:

```
dI/dt = A I + B U        Y = C I
```

- `I` holds the 12 circuit currents and the passive currents.
- `U` holds the 14 power-supply voltages, in the order the controller sends them:
  `CS3U CS2U CS1U CS1L CS2L CS3L PF1 PF2 PF3 PF4 PF5 VS1 PF6 VS3`.
- `Y` holds the plasma outputs.

The plant sends back the 14 coil currents, built from the circuit currents. Every model
therefore behaves as NICE, which has the circuits inside: the controller needs no change
when the full NICE model replaces a linear one.

The actor uses B as it is and discretises the continuous matrices with its time step
(zero-order hold).

## Available models

| Folder | States | Plasma | Vertically unstable | Coordinates |
|---|---|---|---|---|
| `vertical_circuit_105073_t84.3176` | 114 (12 circuits + 102 passive loops) | rigid vertical displacement, Ip fixed | yes | deviation |
| `niceSimpleLin_105073_t84.4176` | 16 (12 circuits + 4 passive) | NICE linearisation, reduced | no | deviation |
| `niceLin_105073_t84.4176` | 1191 (12 circuits + 1179 passive) | NICE linearisation, full | yes | deviation |
| `plasmaless_model` | 114 (12 circuits + 102 passive loops) | no | no | absolute |

The plasma models are all built around pulse 105073 at 84.3–84.4 s (flattop). For other
pulses, build a new model (see below). `plasmaless_model` does not depend on the pulse.

## Using a model

In `cases/overrides/linear_plant_controller_<shot>.ymmsl`:

```yaml
plant.linear_model_dir: "${PDS_REPO}/local_installs/plasmaless-tokamak-circuits/linear_models/vertical_circuit_105073_t84.3176"
plant.operating_point: descriptor   # or f_init
```

- `operating_point` sets where the coil currents and plasma outputs start:
  - `descriptor` starts from the point stored in `linear_model.json`;
  - `f_init` starts from the initial state of the run.
- The power-supply voltages start from the descriptor (zero).
- `plasmaless_model` has absolute coordinates and no operating point:
  - it starts from the coil currents of the initial state, projected on the circuits;
  - it stops with an error if CS1U ≠ CS1L or VS3U ≠ −VS3L;
  - the passive currents start at zero, and `operating_point` is ignored.

## Contract for every model

The actor (`imas_model/configure_linear_model.m`) loads a folder only if:

- `inputs` are the 14 power supplies, in the order above;
- `coil_states` maps the 14 coils to the ITER circuits:
  - CS1U and CS1L on one state with the same sign (series, e.g. both `CS1`);
  - VS3U and VS3L on one state with opposite signs (anti-series, e.g. `VS3` and
    `-VS3`, where `-name` means minus the state `name`);
  - every other coil on a state of its own.

Otherwise it stops with `plasmaless:circuits` or `plasmaless:descriptor`.

`test_all_models_circuits` in `tests/test_plasmaless_linear.m` also checks every folder
here:
- the B columns of CS1U and CS1L are equal;
- VS1 acts as −PF2 −PF3 +PF4 +PF5;
- a positive voltage raises the current of the supply's own coil, and a positive VS3 raises
  the VS3L current;
- after one step, CS1U = CS1L and VS3U = −VS3L.

A new model must pass both: put its folder here and run the tests.

## Folder content

- `linear_model.json` is the descriptor read by the actor. Its main fields are:
  - `files`: the A, B, C matrices (`.npy`);
  - `inputs`: the 14 power supplies, in the order above;
  - `states`, `coil_states`, `passive_states`: which state is which circuit or loop;
  - `outputs`: any of `rmag`, `zmag`, `rX`, `zX`, `Ip`, `gap<k>`, `Dist2ndSep`;
  - `operating_point`: the linearisation point;
  - `convention`: sign changes from the provider's frame to IMAS;
  - `coordinates`: `deviation` (default) or `absolute`.
- `A`, `B`, `C`: continuous-time matrices.

## Raw sources (`sources/`)

These are the models as delivered or exported, before the reduction to circuits. The actor
refuses them.

| Folder | Content | Circuit model |
|---|---|---|
| `sources/niceSimpleLin_105073_t84.4176` | NICE delivery, 17 states (VS3 and -VS3 separate) | `niceSimpleLin_105073_t84.4176` |
| `sources/niceLin_105073_t84.4176` | NICE delivery, 1192 states (VS3 and -VS3 separate) | `niceLin_105073_t84.4176` |
| `sources/plasmaless_model` | 14 independent coils + 102 passive loops, coil-voltage inputs | `plasmaless_model` |

## Building a model (`utils/`)

| Script | Builds |
|---|---|
| `vertical_circuit/make_vertical_circuit_model.m` | the circuit vertical model around an equilibrium (e.g. a NICE inverse run at the case's t0). Set the equilibrium, time and folder name at the top and run it in MATLAB: it writes the folder, reads it back as the plant does and runs the plant tests |
| `make_nicelin_descriptor.py` | `sources/niceLin_…` from the delivery of the NICE full linear model (`linModA/B/C.npy`) |
| `transform_nice_lin_model.m` | `niceLin_…` and `niceSimpleLin_…` from their sources. The states VS3 and -VS3 become one circuit state; this is exact because NICE already holds the anti-series constraint in A and B |
| `export_plasmaless_model.m` | `sources/plasmaless_model` from a machine description |
| `transform_plasmaless_model.m` | `plasmaless_model` from its source, by reduction to the circuits: `Wm = (W'MW)⁻¹W'M`, `A = Wm A0 W`, `B = Wm B0 V` |

`vertical_circuit/` holds everything the circuit model needs:
- a copy of F. Pesamosca's circuit code (plasma-tokamak-circuits, branch
  `vertical_instability`, commit c02e88b);
- `extract_equilibrium_data.m`, the equilibrium reader extended for the PDS entries;
- `write_vertical_circuit_model.m`, the exporter called by the launcher.

Run the scripts on a compute node (MATLAB with IMAS-MATLAB loaded). Each file's header
explains its inputs.
