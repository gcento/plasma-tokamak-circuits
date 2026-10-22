% make_nice_models.m
% Builds the NICE models of linear_models/ from the raw NICE deliveries of
% generators/nice/sources (13 coil states, VS3 and -VS3 separate):
% generators/nice/transform_nice_lin_model.m merges VS3 and -VS3 into one circuit state
% (exact, NICE holds the antiseries constraint in A and B). Then it prints what the plant
% sees in each model (utils/model_summary.m) and, if compare_root is set, checks that it
% reproduces the model of the same name there (utils/same_model.m).
% (generators/nice/make_nicelin_descriptor.py made the raw niceLin source from the
% delivery once; see generators/README.md.)
% Edit the parameters below, then run it in MATLAB on a compute node with an IMAS-MATLAB
% module loaded (e.g. module load PCS/2026.08-pds IMAS-MATLAB/5.6.0-intel-2025b-DD-4.1.1).

%% parameters
here = fileparts(mfilename('fullpath'));
repo_root = fileparts(here);                            % plasmaless-tokamak-circuits
names = {'niceLin_105073_t84.4176', 'niceSimpleLin_105073_t84.4176'};
source_root = fullfile(repo_root,'generators','nice','sources');
out_root = fullfile(repo_root,'linear_models');
overwrite = false;                      % true: replace existing model folders
compare_root = '';                      % '' or a folder holding models of the same names
% Check that the code reproduces the stored models (as tests/test_generators.m):
% out_root = fullfile('/scratch/users',getenv('USER'),'linear_model','regen');
% compare_root = fullfile(repo_root,'linear_models');

%% build, read back, compare
addpath(fullfile(repo_root,'generators','nice'), fullfile(repo_root,'plant'), here);
[~, pf_active_md] = load_machine_description();
for k = 1:numel(names)
    out_dir = fullfile(out_root, names{k});
    transform_nice_lin_model(fullfile(source_root, names{k}), out_dir, overwrite);
    model_summary(out_dir, pf_active_md);
    if ~isempty(compare_root)
        same_model(out_dir, fullfile(compare_root, names{k}));
    end
end
