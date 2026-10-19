function write_state_matrix(out_dir, equilibrium_uri, equilibrium_time, opts)
% Build the vertical-instability model around an equilibrium and write it as a
% continuous-time linear model folder (linear_model.json + A, B, C .npy) for the plant
% setting linear_model_dir, like the other models of linear_models/: the files hold the
% continuous matrices and the plant discretises them with plant.dt (zero-order hold).
%   write_state_matrix(out_dir, equilibrium_uri, equilibrium_time)
%   write_state_matrix(out_dir, equilibrium_uri, equilibrium_time, opts)
% opts (optional fields): passive_currents ('equilibrium', the default, or 'zero' for
% entries without passive currents, e.g. the PDS scenario entries), interpolation
% (ids_get_slice mode, default 3), md_uri (machine description, default data/md_dd4 of
% the plasmaless repository, the plant's default, so that coils and loops match).
% Needs IMAS-MATLAB, the plasmaless repository (add_plasmaless_paths) and this folder on
% the path.
%
% The equilibrium entry needs profiles_2d (psi, j_phi) on a rectangular grid. NICE writes it
% with outputCartGrid = 1, nr_cartgrid, nz_cartgrid in its parameter file (as in
% cases/overrides/inverse_convergence_105073_config_nice_inverse.xml); DINA entries have it.
% Use the same NICE inverse run as preprocess.source_uri of the case, at the case's t0:
% the model's operating point is then the run's initial state.
%
% Example (pulse 105073 at t0 = 84.3176 s of the inverse_convergence run; its pf_passive
% has no passive currents):
%   run(fullfile(getenv('PDS_REPO'),'local_installs','plasmaless-tokamak-circuits','add_plasmaless_paths.m'))
%   addpath(fullfile(getenv('PDS_REPO'),'linear_models','utils','lin_vertical_instability'))
%   write_state_matrix(fullfile(getenv('PDS_REPO'),'linear_models','vertical_105073_t84.3176'), ...
%       fullfile(getenv('PDS_REPO'),'cases','runs','inverse_convergence_105073','out_nice'), ...
%       84.31759781, struct('passive_currents','zero'))
%
% Written model (continuous time, IMAS convention, deviations from the equilibrium):
%   states   the coil currents and the passive loop currents, one state each;
%   inputs   the coil voltages (no input_map);
%   outputs  Ip and rmag: zero rows (the model keeps Ip fixed and has no radial motion);
%            zmag: zIp/Ip0, the rigid vertical displacement of the plasma current, which
%            moves the magnetic axis by the same amount;
%   operating point: coil currents = equilibrium currents Ic0; voltages = R*Ic0, the
%            voltages that keep these currents constant (zero on the coils without
%            resistance in the machine description). No output values: the plant takes Ip
%            and the magnetic axis from the F_INIT equilibrium at t0.
% The time step dt_check below is not in the files: it only discretises the model for the
% read-back check (folder read with configure_linear_model, the plant's loader: same
% discrete A, B and zmag row, same currents and zmag over 1000 steps).

if nargin < 3
    error('write_state_matrix:nargin', ...
        'Use write_state_matrix(out_dir, equilibrium_uri, equilibrium_time[, opts]).');
end
if nargin < 4 || isempty(opts)
    opts = struct();
end
out_dir = char(out_dir);
md_uri = '';
if isfield(opts,'md_uri')
    md_uri = char(opts.md_uri);
    opts = rmfield(opts,'md_uri');
end

% machine description and equilibrium (the extractor matches the same passive loops)
[em_coupling, pf_active, pf_passive] = load_plasmaless_md(md_uri);
opts.pf_passive_md = pf_passive;
equilibrium = extract_equilibrium_data(equilibrium_uri, equilibrium_time, opts);

% continuous model: Ac, Bc and Cc do not depend on the time step; dt_check only sets the
% model's own discrete copy (model.A, model.B), used by the check below
dt_check = 0.002;
model = configure_vertical_instability_model(dt_check, equilibrium, em_coupling, ...
    pf_active, pf_passive);
nc = model.n_coils;
nl = model.n_loops;
n = nc + nl;
coil_names = model.coil_names(:)';
loop_names = model.loop_names(:)';
% configure_linear_model reads coil names as JSON object keys (MATLAB field names)
if ~all(cellfun(@isvarname, coil_names))
    error('write_state_matrix:coilName','Coil names must be valid MATLAB names (got {%s}).', ...
        strjoin(coil_names,','));
end

% outputs: Ip and rmag fixed (zero rows), zmag = zIp/Ip0 (rigid vertical displacement)
y_names = {'Ip','rmag','zmag'};
C = [zeros(2,n); model.Cc(model.ind_zIp,:)/model.Ip0];

% operating point: equilibrium coil currents and the voltages that keep them constant
Ic0 = model.Ic0(:);
V0 = full(model.R(1:nc,:))*Ic0;
if any(Ic0(nc+1:end) ~= 0)
    warning('write_state_matrix:passive', ['The equilibrium has passive currents (max |I| ', ...
        '%.4g A): the descriptor has no passive operating point, so the plant reports only ', ...
        'their deviations.'], max(abs(Ic0(nc+1:end))));
end
coil_states = struct();
current = struct();
voltage = struct();
for ii = 1:nc
    coil_states.(coil_names{ii}) = coil_names{ii};
    current.(coil_names{ii}) = Ic0(ii);
    voltage.(coil_names{ii}) = V0(ii);
end

ev = eig(model.Ac);
unstable = real(ev) > 1e-9;
growth = max(real(ev));
if any(unstable)
    stability_note = sprintf(['Open-loop unstable: %d eigenvalue(s) of A with positive real ', ...
        'part, growth rate %.4g 1/s (e-folding time %.1f ms).'], nnz(unstable), growth, 1e3/growth);
else
    stability_note = 'Open-loop stable: no eigenvalue of A with positive real part.';
end

src = model.source;
shot = 0;
tok = regexp(src.uri,'/(\d{5,6})(/|$)','tokens','once');
if ~isempty(tok)
    shot = str2double(tok{1});
end
d.description = sprintf(['vertical-instability model around %s: rigid vertical plasma ', ...
    'displacement (Ip fixed, no radial motion) with %d coils and %d passive loops; ', ...
    'continuous time, dI/dt = A I + B U, deviations from the equilibrium currents'], ...
    src.description, nc, nl);
d.shot = shot;
d.time = src.time;
d.files = struct('A','A.npy','B','B.npy','C','C.npy');
d.states = [coil_names, loop_names];
d.inputs = coil_names;
d.coil_states = coil_states;
d.passive_states = loop_names;
d.outputs = y_names;
d.coordinates = 'deviation';
d.operating_point = struct('coil_current',current,'voltage',voltage);
d.units = struct('currents','A','voltages','V','rmag_zmag','m','Ip','A');
d.provenance = sprintf(['configure_vertical_instability_model (A. Coutlis et al., Nucl. ', ...
    'Fusion 39 (1999); F. Pesamosca, PhD thesis, EPFL (2021), App. B); passive currents: %s; ', ...
    'written by linear_models/utils/lin_vertical_instability/write_state_matrix.m on %s'], ...
    src.passive_source, char(datetime('today','Format','yyyy-MM-dd')));
d.notes = { ...
    ['Continuous-time matrices, as for the other models: the plant discretises them with ', ...
     'plant.dt (zero-order hold).'], ...
    ['zmag = zIp/Ip0: rigid vertical displacement of the plasma current distribution; the ', ...
     'magnetic axis moves by the same amount. The Ip and rmag rows are zero: Ip is fixed and ', ...
     'there is no radial motion, so the controller cannot act on them.'], ...
    ['Operating point: coil currents = equilibrium currents Ic0; voltages = R*Ic0, which keep ', ...
     'the currents constant (zero on the coils without resistance in the machine ', ...
     'description). No output values: with operating_point descriptor the plant takes Ip and ', ...
     'the magnetic axis from the F_INIT equilibrium at t0.'], ...
    stability_note};

if ~isfolder(out_dir)
    mkdir(out_dir);
end
write_npy(fullfile(out_dir,'A.npy'), model.Ac);
write_npy(fullfile(out_dir,'B.npy'), model.Bc);
write_npy(fullfile(out_dir,'C.npy'), C);
fid = fopen(fullfile(out_dir,'linear_model.json'),'w');
if fid < 0
    error('write_state_matrix:file','Cannot write %s.', fullfile(out_dir,'linear_model.json'));
end
fprintf(fid,'%s\n',jsonencode(d,'PrettyPrint',true));
fclose(fid);

% check: the folder read by the plant's own loader, discretised with dt_check, against the
% model's own discrete copy
w = warning('off','plasmaless:unstable');   % expected for this model, reported below
restore_warning = onCleanup(@() warning(w));
lin = configure_linear_model(out_dir, dt_check, pf_active);
clear restore_warning

rel = @(a,b) max(abs(a(:)-b(:)))/max(abs(b(:)));
z_row = C(3,:);
err_matrices = max([rel(lin.A, model.A), rel(lin.B, model.B), rel(lin.C_y(3,:), z_row)]);
if err_matrices > 1e-10
    error('write_state_matrix:check', ['the folder read back differs from the model: max ', ...
        'relative difference of A, B and the zmag row %g.'], err_matrices);
end
rng_state = rng;
rng(0);
x_model = zeros(n,1);
x_lin = zeros(n,1);
diff_I = 0; max_I = 0; diff_z = 0; max_z = 0;
for k = 1:1000
    u = 100*randn(nc,1);
    x_model = model.A*x_model + model.B*u;
    x_lin = lin.A*x_lin + lin.B*u;
    I_lin = [lin.x_to_coil*x_lin; lin.x_to_loop*x_lin];
    diff_I = max(diff_I, max(abs(I_lin - x_model)));
    max_I = max(max_I, max(abs(x_model)));
    diff_z = max(diff_z, abs(lin.C_y(3,:)*x_lin - z_row*x_model));
    max_z = max(max_z, abs(z_row*x_model));
end
rng(rng_state);
err_run = max(diff_I/max_I, diff_z/max_z);
if err_run > 1e-9
    error('write_state_matrix:check', ['currents or zmag of the folder read back differ from ', ...
        'the model over 1000 steps (max relative difference %g).'], err_run);
end

fprintf(['write_state_matrix: wrote %s\n  continuous A, B, C: %d states (%d coils + %d ', ...
    'passive loops), outputs Ip, rmag, zmag; the plant discretises them with plant.dt\n', ...
    '  check (dt %g): discrete A, B and zmag row as the model (max relative difference ', ...
    '%.1e); currents and zmag equal over 1000 steps (%.1e)\n  %s\n'], ...
    out_dir, n, nc, nl, dt_check, err_matrices, err_run, stability_note);
nz = find(V0 ~= 0)';
if ~isempty(nz)
    parts = arrayfun(@(k) sprintf('%s %.3g V', coil_names{k}, V0(k)), nz, 'UniformOutput', false);
    fprintf('  operating-point voltages R*Ic0, non-zero entries: %s\n', strjoin(parts, ', '));
end
fprintf('  use it with: plant.linear_model_dir: "%s"\n', out_dir);
end

function write_npy(filename, M)
% minimal .npy writer: float64, 2-D, Fortran order (MATLAB's column-major storage)
header = sprintf('{''descr'': ''<f8'', ''fortran_order'': True, ''shape'': (%d, %d), }', ...
    size(M,1), size(M,2));
pad = mod(-(10 + numel(header) + 1), 64);    % total header length a multiple of 64
header = [header repmat(' ',1,pad) newline];
fid = fopen(filename,'w','ieee-le');
fwrite(fid,[147 double('NUMPY') 1 0],'uint8');
fwrite(fid,numel(header),'uint16');
fwrite(fid,header,'char');
fwrite(fid,full(M(:)),'float64');
fclose(fid);
end
