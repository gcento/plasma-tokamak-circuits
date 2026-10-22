# Model generators

How the models of `../linear_models/` are made. Every generator writes a model folder
(`linear_model.json` + continuous A, B, C) with `common/write_linear_model.m`, which reads
the matrix files back and refuses to overwrite an existing model unless asked. The
launchers with the parameters at the top are in `../utils/`. Run them in MATLAB on a
compute node with the IMAS-MATLAB module given in their header.

| Model | From | Code | Launcher |
|---|---|---|---|
| `plasmaless_model` | machine description | `plasmaless/` | `utils/make_plasmaless_model.m` |
| `vertical_circuit_<shot>_t<time>` | machine description + equilibrium | `vertical/` | `utils/make_vertical_circuit_model.m` |
| `niceLin_…`, `niceSimpleLin_…` | NICE deliveries | `nice/` | `utils/make_nice_models.m` |
| (`../data/md_dd4`, the machine description) | DD3 IDSs | `md/convert_md_dd3_to_dd4.py` | `python md/convert_md_dd3_to_dd4.py` (IMAS-MUSCLE3 module) |

## Plasma-free model (`plasmaless/`)

Two steps:

1. `export_plasmaless_model.m`: the 14 coils and the 102 passive loops of the machine
   description as independent circuits, with the 14 coil voltages as inputs:
   A0 = −M⁻¹R, B0 = M⁻¹[I14; 0] (M, R: single-coil mutual inductances and resistances
   of the machine description). Output: the raw source `plasmaless/source/`, which the
   plant refuses because its inputs are coil voltages.
2. `transform_plasmaless_model.m`: Galerkin reduction to the 12 ITER circuits.
   - x = W x_c with W = blkdiag(Wa, I): CS1U = CS1L = CS1, VS3U = VS3, VS3L = −VS3.
   - Coil voltages from the 14 supplies: u = V u_ps (VS1 on PF2..PF5 with −, −, +, +;
     VS3U = −Vps/2, VS3L = +Vps/2).
   - Wm = (W'MW)⁻¹W'M, A = Wm A0 W, B = Wm B0 V, C = C0 W. W' sums the coil equations
     along each series loop and cancels the constraint voltages.

   Output: `linear_models/plasmaless_model` (absolute coordinates).

The result is A = −(W'MW)⁻¹W'RW and B = (W'MW)⁻¹[Wa'V; 0], with Wa'V the supply-to-circuit
map. These are the circuit matrices that `common/configure_plasmaless_model.m` builds
for the vertical model (M11 there), up to the sign of the VS3 state (VS3L = VS3 there).

## Vertical model (`vertical/`)

Rigid vertical displacement of the plasma current around an equilibrium, Ip fixed, no
radial motion (A. Coutlis et al., Nucl. Fusion 39 (1999); F. Pesamosca, PhD thesis, EPFL
(2021), App. B):

1. `extract_equilibrium_data.m` reads the equilibrium slice:
   - (R, Z) grid, toroidal current density, Ip0, closed LCFS;
   - the coil currents;
   - psi sign checked with Ampère's law;
   - passive currents from the entry or zero (the PDS entries have none).
2. `configure_vertical_instability_model.m` builds the circuit model:
   - circuit inductance M11 and resistance from the machine description
     (`common/configure_plasmaless_model.m`, circuit version);
   - plasma–conductor couplings M12, M21 from the vertical derivative of the current
     density and the radial field of the conductors;
   - vertical force balance M22, eliminated with a Schur complement:
     Meff = M11 − M12 M21/M22, A = −Meff⁻¹R, B = Meff⁻¹ S P.
3. `write_vertical_circuit_model.m` writes the folder.
   - Outputs: Ip and rmag (zero rows), zmag = zIp/Ip0.
   - Operating point: the equilibrium coil currents.
   - It checks the input order and the VS1 column (−PF2 −PF3 +PF4 +PF5).

Launcher parameters: the equilibrium (e.g. the `out_nice` of an `inverse_convergence` run
of IMAS-PDS), the time, the model name and the passive currents.

## NICE models (`nice/`)

`nice/sources/<name>` holds the delivered linear models (provider's frame, 13 coil
states with VS3 and −VS3 separate, A(:,−VS3) = −A(:,VS3)) with their descriptors; the
README of `niceSimpleLin_105073_t84.4176` records the delivery, the pulse correction and
the derivation of the sign frame.

- `make_nicelin_descriptor.py` made the raw niceLin source once from the delivery of the
  full model (C rows rX2, zX2 dropped; descriptor copied from niceSimpleLin). It needs the
  original delivery folder, not in this repository.
- `transform_nice_lin_model.m` merges VS3 and −VS3 into one circuit state:
  VS3 = (VS3 − (−VS3))/2, A = T A0 W, B = T B0, C = C0 W. This is exact because NICE holds
  the anti-series constraint in A and B.

## Shared code (`common/`)

- From F. Pesamosca's plasma-tokamak-circuits (branch `vertical_instability`, commit
  c02e88b):
  - `configure_plasmaless_model.m`: circuit version, the 12 circuits + passive loops with
    the 14 supplies as inputs; also gives the single-coil M and R;
  - `plasmaless_circuit_topology.m`: Ia = T Ic, Vc = P Vps;
  - `map_coils_to_circuits.m`, `map_powersupplies_to_circuits.m`.

  Renaming the plant code from "plasmaless" to "linear_model" (2026-10-22) changed only
  one name in this copy and in `vertical/`: the machine-description loader, now
  `plant/load_machine_description.m` (formerly `load_plasmaless_md.m`).
- IMAS readers of the machine description: `get_coupling_matrix.m`,
  `plasmaless_grid_from_ids.m`.
- Writers: `write_npy.m` (read back by `plant/read_npy.m`), `write_linear_model.m`.

## Checking that the stored models are reproduced

- `tests/test_generators.m` regenerates the plasma-free source and model and the two NICE
  models in `/scratch/users/$USER/linear_model` and compares them with the stored files
  (matrices within 1e-12, descriptor fields equal).
- Vertical model: in `utils/make_vertical_circuit_model.m`, set the run the stored model
  was built from (`inverse_convergence_105073_20261017_121544`), `out_dir` in scratch and
  `compare_with` to the stored folder. `utils/same_model.m` then makes the same check as
  the test and prints one line.
