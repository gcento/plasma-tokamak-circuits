% make_vertical_circuit_model.m
% Builds the vertical-instability model with circuit states (rigid vertical plasma
% displacement around an equilibrium, open-loop unstable) and writes it as a model folder
% for the PDS plant (linear_plant_controller): generators/vertical/
% write_vertical_circuit_model.m (generators/README.md). Then it prints what the plant
% sees in the folder (utils/model_summary.m), checks that it reproduces the folder
% compare_with if set (utils/same_model.m) and runs the plant tests.
% Edit the parameters below, then run it in MATLAB on a compute node with an IMAS-MATLAB
% module loaded (e.g. module load PCS/2026.08-pds IMAS-MATLAB/5.6.0-intel-2025b-DD-4.1.1).

%% parameters
here = fileparts(mfilename('fullpath'));
repo_root = fileparts(here);                            % plasmaless-tokamak-circuits
pds_repo = fileparts(fileparts(repo_root));             % IMAS-PDS (local_installs/..)

% equilibrium with profiles_2d (IMAS URI or entry folder): the out_nice of the
% inverse_convergence run, at the case's t0. A folder path is resolved (cases/runs/<case>
% is a symlink to the latest run), so the descriptor records the run actually used.
equilibrium = fullfile(pds_repo,'cases','runs','inverse_convergence_105073','out_nice');
t0 = 84.31759781;                       % [s]
model_name = 'vertical_circuit_105073_t84.3176';
passive_currents = 'zero';              % the PDS entries have no passive currents
out_dir = fullfile(repo_root,'linear_models',model_name);
overwrite = false;                      % true: replace an existing out_dir
compare_with = '';                      % '' or a model folder that out_dir must reproduce
run_tests = true;                       % tests/test_linear_models.m at the end
% Check that the stored model is reproduced (it was built on 2026-10-20 from the run
% inverse_convergence_105073_20261017_121544):
% equilibrium = fullfile(pds_repo,'cases','runs','inverse_convergence_105073_20261017_121544','out_nice');
% out_dir = fullfile('/scratch/users',getenv('USER'),'linear_model','regen',model_name);
% compare_with = fullfile(repo_root,'linear_models',model_name);

%% 1) build and write the model folder
if isfolder(equilibrium)
    equilibrium = char(java.io.File(equilibrium).getCanonicalPath());
end
addpath(fullfile(repo_root,'generators','vertical'));
write_vertical_circuit_model(out_dir, equilibrium, t0, ...
    struct('passive_currents',passive_currents,'overwrite',overwrite));

%% 2) read it back as the PDS plant does, compare
addpath(fullfile(repo_root,'plant'), here);
[~, pf_active_md] = load_machine_description();
model_summary(out_dir, pf_active_md);
if ~isempty(compare_with)
    same_model(out_dir, compare_with);
end

%% 3) plant tests (they also load every folder of linear_models/)
if run_tests
    results = runtests(fullfile(repo_root,'tests','test_linear_models.m'));
    disp(table(results));
end
