% make_vertical_circuit_model.m
% Writes the circuit vertical model folder used by the PDS plant (linear_plant_controller),
% reads it back with the plant's own loader and runs the plant tests.
% Edit the parameters below, then run this script in MATLAB on a compute node with the
% IMAS-MATLAB module loaded (e.g. module load PCS/2026.08-pds
% IMAS-MATLAB/5.6.0-intel-2025b-DD-4.1.1).

%% parameters
here = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(here)));      % plasmaless-tokamak-circuits
pds_repo = fileparts(fileparts(repo_root));             % IMAS-PDS (local_installs/..)

% equilibrium with profiles_2d: the out_nice of the inverse_convergence run, at the
% case's t0
equilibrium = fullfile(pds_repo,'cases','runs','inverse_convergence_105073','out_nice');
t0 = 84.31759781;                       % [s]
model_name = 'vertical_circuit_105073_t84.3176';
passive_currents = 'zero';              % the PDS entries have no passive currents
run_tests = true;                       % tests/test_plasmaless_linear.m at the end

%% 1) build and write the model folder
out_dir = fullfile(repo_root,'linear_models',model_name);
addpath(here);
write_vertical_circuit_model(out_dir, equilibrium, t0, ...
    struct('passive_currents',passive_currents));
rmpath(here);   % its configure_plasmaless_model.m (circuit version) must not shadow the plant's

%% 2) read it back as the PDS plant does
addpath(repo_root, fullfile(repo_root,'imas_model'));
[~, pf_active_md] = load_plasmaless_md();
w = warning('off','plasmaless:unstable');       % open-loop unstable by design
model = configure_linear_model(out_dir, 1e-3, pf_active_md);
warning(w);
row = @(name) model.x_to_coil(strcmp(model.coil_names,name),:);
fprintf(['plant loader: %d states (%d coils, %d passive loops), %d inputs {%s}\n', ...
    '  CS1U = CS1L: %d, VS3U = -VS3L: %d, growth rate %.4g 1/s\n'], model.n_states, ...
    model.n_coils, model.n_loops, numel(model.input_names), strjoin(model.input_names,' '), ...
    isequal(row('CS1U'),row('CS1L')), isequal(row('VS3U'),-row('VS3L')), ...
    max(real(eig(model.A_c))));

%% 3) plant tests (they also load this model folder)
if run_tests
    results = runtests(fullfile(repo_root,'tests','test_plasmaless_linear.m'));
    disp(table(results));
end
