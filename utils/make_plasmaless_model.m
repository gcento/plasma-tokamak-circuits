% make_plasmaless_model.m
% Builds the plasma-free circuit model linear_models/plasmaless_model from the machine
% description, in two steps (generators/plasmaless, generators/README.md):
%   1) export_plasmaless_model: the raw coil model, 14 independent coils + passive loops
%      with the coil voltages as inputs, A0 = -M\R, B0 = M\[I14; 0]
%      (generators/plasmaless/source; optional, run_export);
%   2) transform_plasmaless_model: the Galerkin reduction to the 12 ITER circuits with the
%      14 power-supply voltages as inputs, Wm = (W'MW)\(W'M), A = Wm A0 W, B = Wm B0 V.
% Then it prints what the plant sees in the model (utils/model_summary.m) and, if
% compare_with is set, checks that the model reproduces that folder (utils/same_model.m).
% Edit the parameters below, then run it in MATLAB on a compute node with an IMAS-MATLAB
% module loaded (e.g. module load PCS/2026.08-pds IMAS-MATLAB/5.6.0-intel-2025b-DD-4.1.1).

%% parameters
here = fileparts(mfilename('fullpath'));
repo_root = fileparts(here);                            % plasmaless-tokamak-circuits
md_uri = '';                            % machine description ('' = data/md_dd4)
run_export = false;                     % true: rebuild the raw source first (step 1)
source_dir = fullfile(repo_root,'generators','plasmaless','source');
out_dir = fullfile(repo_root,'linear_models','plasmaless_model');
overwrite = false;                      % true: replace existing source_dir/out_dir models
compare_with = '';                      % '' or a model folder that out_dir must reproduce
% Check that the code reproduces the stored files (the same check as tests/test_generators.m):
% scratch = fullfile('/scratch/users',getenv('USER'),'linear_model','regen');
% run_export = true; source_dir = fullfile(scratch,'plasmaless_source');
% out_dir = fullfile(scratch,'plasmaless_model');
% compare_with = fullfile(repo_root,'linear_models','plasmaless_model');

%% 1) raw coil model, 2) circuit model
addpath(fullfile(repo_root,'generators','plasmaless'), fullfile(repo_root,'plant'), here);
if run_export
    export_plasmaless_model(source_dir, md_uri, overwrite);
end
transform_plasmaless_model(source_dir, out_dir, md_uri, overwrite);

%% 3) read it back as the PDS plant does, compare
[~, pf_active_md] = load_machine_description(md_uri);
model_summary(out_dir, pf_active_md);
if ~isempty(compare_with)
    same_model(out_dir, compare_with);
end
