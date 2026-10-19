# Linear plant models

Linear models of the ITER coils, vessel and (optionally) plasma. In IMAS-PDS they are selected by the workflow `linear_plant_controller`.

Each folder is one model:

```
dI/dt = A I + B U        Y = C I
```

- `I` holds the coil and passive currents.
- `U` holds the coil voltages.
- `Y` holds the plasma outputs.

The actor discretises the continuous matrices with its time step (zero-order hold).

## Available models

| Folder | States | Plasma | Vertically unstable | |
|---|---|---|---|---|
| `plasmaless_model` | 116 (14 coils + 102 passive loops) | no | no |
| `vertical_unstable_model_105073_t84.3176` | 116 | rigid vertical displacement, Ip fixed | yes |
| `niceSimpleLin_105073_t84.4176` | 17 | NICE linearisation, reduced | no |
| `niceLin_105073_t84.4176` | 1192 | NICE linearisation, full | yes |

The plasma models are all built around pulse 105073 at about 84.4 s (flattop). For other
pulses, build a new model (see below) or use `plasmaless_model`.


## Using a model

In `cases/overrides/linear_plant_controller_<shot>.ymmsl`:

```yaml
plant.linear_model_dir: "${PDS_REPO}/local_installs/plasmaless-tokamak-circuits/linear_models/niceLin_105073_t84.4176"
plant.operating_point: descriptor   # or f_init
```

- An empty `linear_model_dir` makes the actor build the plasmaless circuit model itself.
- `operating_point` matters only for the plasma models, which work in deviations from the
  linearisation point:
  - `descriptor` starts from the point stored in `linear_model.json`;
  - `f_init` starts from the initial state of the run, coil voltages included.
- `plasmaless_model` uses absolute currents and always starts from the initial coil
  currents.

## Folder content

- `linear_model.json` is the descriptor read by the actor
  (`imas_model/configure_linear_model.m`). Its main fields are:
  - `files`: the A, B, C matrices (`.npy`);
  - `states`, `coil_states`, `passive_states`: which state is which coil or loop;
  - `inputs`, `input_map`: how each model input is made of the 14 coil voltages;
  - `outputs`: any of `rmag`, `zmag`, `rX`, `zX`, `Ip`, `gap<k>`, `Dist2ndSep`;
  - `coordinates`: `deviation` (default) or `absolute`;
  - `operating_point`: the linearisation point (deviation models only);
  - `convention`: sign changes from the provider's frame to IMAS.
- `A`, `B`, `C`: continuous-time matrices.

The actor takes the 14 voltages in coil order (CS3U … PF6, VS3U, VS3L). Generator
voltages are expressed through `input_map`:
- VS3 is the anti-series pair, `VS3 = V(VS3L) − V(VS3U)`.
- VS1 has no coil of its own. The controller adds it to PF2..PF5 with signs − − + +,
  which matches the VS1 column of B in the NICE models.

## Building a new model (`utils/`)

| Script | Builds |
|---|---|
| `lin_vertical_instability/write_state_matrix.m` | a vertical-instability model around any equilibrium (e.g. a NICE inverse run at the case's t0) |
| `export_plasmaless_model.m` | `plasmaless_model` from a machine description |
| `make_nicelin_descriptor.py` | the folder of a NICE full linear model (`linModA/B/C.npy`) |
