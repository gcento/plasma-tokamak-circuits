# Linear plant models

Linear models of the ITER coils, vessel and plasma. In IMAS-PDS they are selected by the
workflow `linear_plant_controller`.

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
| `vertical_circuit_105073_t84.3176` | 114 (12 circuits + 102 passive loops) | rigid vertical displacement, Ip fixed | yes (7.47 1/s) | deviation |
| `niceSimpleLin_105073_t84.4176` | 16 (12 circuits + 4 passive) | NICE linearisation, reduced | no | deviation |
| `niceLin_105073_t84.4176` | 1191 (12 circuits + 1179 passive) | NICE linearisation, full | yes | deviation |
| `plasmaless_model` | 114 (12 circuits + 102 passive loops) | no | no | absolute |

The plasma models are all built around pulse 105073 at 84.3–84.4 s (flattop). For other
pulses, build a new model (`../generators/README.md`). `plasmaless_model` does not depend
on the pulse.

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

To see how a model responds to the power supplies before using it in a closed loop:
`../utils/run_supply_sequence.m` (models overlaid, `../README.md`).

## Contract for every model

The actor (`../plant/configure_linear_model.m`) loads a folder only if:

- `inputs` are the 14 power supplies, in the order above;
- `coil_states` maps the 14 coils to the ITER circuits:
  - CS1U and CS1L on one state with the same sign (series, e.g. both `CS1`);
  - VS3U and VS3L on one state with opposite signs (anti-series, e.g. `VS3` and
    `-VS3`, where `-name` means minus the state `name`);
  - every other coil on a state of its own.

Otherwise it stops with `linear_model:circuits` or `linear_model:descriptor`.

`test_all_models_circuits` in `../tests/test_linear_models.m` also checks every folder
here:
- the B columns of CS1U and CS1L are equal;
- VS1 acts as −PF2 −PF3 +PF4 +PF5;
- a positive voltage raises the current of the supply's own coil, and a positive VS3 raises
  the VS3L current;
- after one step, CS1U = CS1L and VS3U = −VS3L.

A new model must pass both: put its folder here and run the tests.

## Folder content

- `A`, `B`, `C`: continuous-time matrices (`.npy`, or whitespace-separated `.txt`).
- `linear_model.json`: the descriptor read by the actor.

| Field | Content |
|---|---|
| `description`, `shot`, `time` | free text and the linearisation point |
| `files` | `A`, `B`, `C`: file names of the matrices |
| `states` | the n state names (rows of A) |
| `inputs` | the 14 power supplies, in the order above (columns of B, used as they are; `input_map` is refused) |
| `coil_states` | coil short name → state name; `-name` (when no state is called `-name`) means minus the state `name` |
| `passive_states` | the passive states; every state is a coil state or a passive state, not both |
| `outputs` | the p output names (rows of C), among `rmag`, `zmag`, `rX`, `zX`, `Ip`, `gap<k>`, `Dist2ndSep` |
| `operating_point` (optional) | `coil_current` (coil → A), `voltage` (supply → V), `outputs` (output → value), absolute values at the linearisation point, IMAS convention |
| `convention` (optional) | sign frame of the matrix files (below) |
| `coordinates` (optional) | `deviation` (default) or `absolute` (no `operating_point` allowed) |
| `notes`, `units`, `provenance`, `provider_note` | free text |

**Convention.** `state_sign`, `input_sign` (±1) and `output_sign` (output → ±1) give the
frame of the matrix files with respect to IMAS: `x_file = state_sign·x`,
`u_file = input_sign·u`, `y_file(k) = output_sign(k)·y(k)`. The loader applies, before the
discretisation, `B = (input_sign/state_sign)·B_file` and
`C(k,:) = (state_sign/output_sign(k))·C_file(k,:)`; A is unchanged. Absent entries are +1.
The NICE models have `state_sign = input_sign = −1` and `output_sign(Ip) = −1`.

**Outputs in the equilibrium IDS** (`../plant/plant_observables_to_equilibrium.m`,
DD 4.1): `rmag`, `zmag` → `global_quantities.magnetic_axis`; `rX`, `zX` →
`contour_tree.node{1}` with `critical_type` 1; `Ip` → `global_quantities.ip`; `gap<k>` →
`boundary.gap{k}.value`; `Dist2ndSep` has no DD4 field and is not sent.

**Errors** of the loader: `linear_model:descriptor` (fields, names, inputs, operating point,
convention), `linear_model:circuits` (coil states not the ITER circuits), `linear_model:size`
(matrix sizes), `linear_model:observable` (output name), `linear_model:file` (missing matrix).
`linear_model:unstable` is a warning when A has an eigenvalue with real part > 1e-9.

## How the models are made

The raw sources (NICE deliveries, the plasma-free coil model) and the generators are in
`../generators/` (`../generators/README.md`); the actor refuses the raw sources.
