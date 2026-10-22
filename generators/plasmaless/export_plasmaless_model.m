function export_plasmaless_model(out_dir, md_uri, overwrite)
% Export the plasmaless coil+vessel model (no plasma) of a machine description as the raw
% source generators/plasmaless/source (linear_model.json + A, B, C .npy): 14 independent
% coils + passive loops, with the 14 coil voltages as inputs. The plant does not load it
% (its inputs are not the power supplies): transform_plasmaless_model.m reduces it to the
% ITER circuits and writes linear_models/plasmaless_model.
%   export_plasmaless_model(out_dir)                     % md: data/md_dd4 of this repository
%   export_plasmaless_model(out_dir, md_uri)
%   export_plasmaless_model(out_dir, md_uri, overwrite)  % overwrite = true replaces out_dir
% Launcher: utils/make_plasmaless_model.m. Needs IMAS-MATLAB (reads the machine
% description) and the Control System Toolbox; run it on a compute node.
%
% M and R are the single-coil matrices of the machine description (M_coils, R_coils of
% generators/common/configure_plasmaless_model.m): A0 = -M\R, B0 = M\[I14; 0]. The model is
% absolute (states = coil and passive currents), so the descriptor has coordinates
% 'absolute' and no operating point.

if nargin < 2
    md_uri = '';
end
if nargin < 3 || isempty(overwrite)
    overwrite = false;
end
here = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(here));                 % generators/plasmaless
addpath(fullfile(repo_root,'plant'), fullfile(repo_root,'generators','common'));
if isempty(md_uri)
    md_uri = ['imas:hdf5?path=' fullfile(repo_root,'data','md_dd4')];
end
md_uri = char(md_uri);   % imas_open needs a char URI

% dt is not in the exported matrices: configure_plasmaless_model needs one (M and R do not
% depend on it)
dt = 0.002;
[em_coupling, pf_active, pf_passive] = load_machine_description(md_uri);
base = configure_plasmaless_model(em_coupling, pf_active, pf_passive, dt);
n_coils = base.n_coils;
n = n_coils + base.n_loops;
M = base.M_coils;
R = base.R_coils;
% continuous matrices of the 14 independent coils and the passive loops
Ac = -M\R;
Bc = M\eye(n,n_coils);

d.description = ['plasmaless coil+vessel circuit model (no plasma), ', num2str(n_coils), ...
    ' coils + ', num2str(base.n_loops), ' passive loops: dI/dt = A I + B U, ', ...
    'absolute currents and voltages'];
d.shot = 0;
d.time = 0;
d.files = struct('A','A.npy','B','B.npy','C','C.npy');
d.states = [base.coil_names base.loop_names];
d.inputs = base.coil_names;
d.coil_states = struct();
for ii = 1:n_coils
    d.coil_states.(base.coil_names{ii}) = base.coil_names{ii};
end
d.passive_states = base.loop_names;
d.outputs = {};
d.coordinates = 'absolute';
d.units = struct('currents','A','voltages','V');
d.provenance = sprintf(['machine description %s; A = -M\\R, B = M\\eye(n,n_coils) with the ', ...
    'single-coil M, R of generators/common/configure_plasmaless_model.m; written by ', ...
    'generators/plasmaless/export_plasmaless_model.m on %s. Raw source: the inputs are the ', ...
    '14 coil voltages, so the plant does not load it; ', ...
    'generators/plasmaless/transform_plasmaless_model.m writes its circuit form ', ...
    'linear_models/plasmaless_model'], md_uri, char(datetime('today','Format','yyyy-MM-dd')));
d.notes = {['Not tied to a pulse: shot and time are 0. Absolute model (no operating point): ', ...
    'the actor starts it from the F_INIT coil currents at t0 (passive currents 0) and ', ...
    'ignores the setting operating_point.']};
write_linear_model(out_dir, Ac, Bc, zeros(0,n), d, overwrite);
fprintf(['export_plasmaless_model: wrote %s (%d states, %d coil-voltage inputs); next: ', ...
    'transform_plasmaless_model.m\n'], out_dir, n, n_coils);
end
