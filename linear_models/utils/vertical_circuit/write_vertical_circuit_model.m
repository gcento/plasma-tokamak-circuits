function write_vertical_circuit_model(out_dir, equilibrium_uri, equilibrium_time, opts)
% Build the vertical-instability model with circuit states (configure_vertical_instability_model:
% 12 active circuits + passive loops, the 14 power-supply voltages as inputs) around an
% equilibrium, and write it as a linear model folder (linear_model.json + continuous A, B,
% C .npy) for the plant setting linear_model_dir. The plant discretises the matrices with
% plant.dt (zero-order hold). Usually run through make_vertical_circuit_model.m.
%   write_vertical_circuit_model(out_dir, equilibrium_uri, equilibrium_time)
%   write_vertical_circuit_model(out_dir, equilibrium_uri, equilibrium_time, opts)
% opts (optional fields):
%   md_uri            machine description (default: data/md_dd4 of this repository, the
%                     plant's default)
%   passive_currents  'equilibrium' (default) or 'zero' (the PDS entries have no passive
%                     currents)
%   interpolation     ids_get_slice mode (default 3)
%   shot              pulse number for the descriptor (default: the first 6-digit number
%                     of equilibrium_uri)
% The model code of this folder is a copy of F. Pesamosca's plasma-tokamak-circuits
% (branch vertical_instability, commit c02e88b): configure_vertical_instability_model.m,
% configure_plasmaless_model.m (circuit version: it shadows the coil version of
% imas_model/ while this function runs), plasmaless_circuit_topology.m,
% map_coils_to_circuits.m and map_powersupplies_to_circuits.m (not called here: the
% supply -> circuit map, Vc = P*Vps); extract_equilibrium_data.m is the extended reader (DD4 psi sign,
% closed LCFS for the PDS entries, zero passive currents). The other helpers come from
% imas_model/. The MATLAB path is restored on exit. Needs IMAS-MATLAB and the Control
% System Toolbox: run it on a compute node.
%
% Written model (continuous time, IMAS convention, deviations from the equilibrium):
%   states   12 circuit currents {CS3U CS2U CS1 CS2L CS3L PF1..PF6 VS3}, then passive loops
%   inputs   14 power-supply voltages {CS3U CS2U CS1U CS1L CS2L CS3L PF1..PF5 VS1 PF6 VS3}
%   coils    CS1U = CS1L = CS1 (series), VS3U = -VS3, VS3L = VS3 (antiseries)
%   outputs  Ip and rmag: zero rows (Ip fixed, no radial motion); zmag = zIp/Ip0
%   operating point: coil currents = equilibrium coil currents; no voltages (zero)

if nargin < 3
    error('write_vertical_circuit_model:nargin', ['Use write_vertical_circuit_model(out_dir, ', ...
        'equilibrium_uri, equilibrium_time[, opts]).']);
end
if nargin < 4 || isempty(opts)
    opts = struct();
end
out_dir = char(out_dir);
[md_uri, opts] = take_option(opts, 'md_uri', '');
[shot, opts] = take_option(opts, 'shot', []);

% paths: imas_model/ of this repository for the shared helpers, this folder on top of it
% (its configure_plasmaless_model.m is the circuit version); restored on exit
here = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(here)));    % linear_models/utils/vertical_circuit
saved_path = path;
restore_path = onCleanup(@() path(saved_path));
addpath(fullfile(repo_root,'imas_model'));
addpath(here);
expected = {
    'configure_vertical_instability_model', here
    'configure_plasmaless_model',           here
    'plasmaless_circuit_topology',          here
    'map_coils_to_circuits',                here
    'extract_equilibrium_data',             here
    'load_plasmaless_md',                   fullfile(repo_root,'imas_model')
    'plasmaless_grid_from_ids',             fullfile(repo_root,'imas_model')};
canon = @(p) char(java.io.File(p).getCanonicalPath());
for k = 1:size(expected,1)
    found = which(expected{k,1});
    if isempty(found) || ~strcmp(canon(fileparts(found)), canon(expected{k,2}))
        error('write_vertical_circuit_model:path', '%s resolves to ''%s'', expected in %s.', ...
            expected{k,1}, found, expected{k,2});
    end
end

% machine description, equilibrium and model (Ac, Bc, Cc do not depend on dt_check)
[em_coupling, pf_active, pf_passive] = load_plasmaless_md(md_uri);
opts.pf_passive_md = pf_passive;
equilibrium = extract_equilibrium_data(equilibrium_uri, equilibrium_time, opts);
dt_check = 0.002;
model = configure_vertical_instability_model(dt_check, equilibrium, em_coupling, ...
    pf_active, pf_passive);

% checks: inputs in the controller's order, VS1 sign of NICE and PCSSP, coil-circuit map
supply_names = {'CS3U','CS2U','CS1U','CS1L','CS2L','CS3L','PF1','PF2','PF3','PF4', ...
    'PF5','VS1','PF6','VS3'};
if ~isequal(model.ps_names(:)', supply_names)
    error('write_vertical_circuit_model:inputs', 'Power supplies {%s} differ from {%s}.', ...
        strjoin(model.ps_names,','), strjoin(supply_names,','));
end
col = @(name) model.Bc(:,strcmp(supply_names,name));
vs1_ref = -col('PF2') - col('PF3') + col('PF4') + col('PF5');
if norm(col('VS1') - vs1_ref) > 1e-10*norm(vs1_ref)
    error('write_vertical_circuit_model:vs1', ['The VS1 column of B is not -PF2 -PF3 +PF4 ', ...
        '+PF5 (the convention of NICE and of the controller): the plasmaless_circuit_topology.m ', ...
        'of %s is older than plasma-tokamak-circuits commit c02e88b.'], here);
end
nc = model.n_circuits;
nl = model.n_loops;
n = nc + nl;
coil_names = model.coil_names(:)';
circuit_names = model.circuit_names(:)';
loop_names = model.loop_names(:)';
if ~all(cellfun(@isvarname, coil_names))
    error('write_vertical_circuit_model:coilName', ['Coil names must be valid MATLAB ', ...
        'names (got {%s}).'], strjoin(coil_names,','));
end
coil_states = struct();
for ii = 1:numel(coil_names)
    k = find(model.T(ii,:));
    if numel(k) ~= 1 || abs(model.T(ii,k)) ~= 1
        error('write_vertical_circuit_model:topology', 'Coil %s is not one circuit with +-1.', ...
            coil_names{ii});
    end
    if model.T(ii,k) > 0
        coil_states.(coil_names{ii}) = circuit_names{k};
    else
        coil_states.(coil_names{ii}) = ['-' circuit_names{k}];
    end
end

% outputs: Ip and rmag fixed (zero rows), zmag = zIp/Ip0
y_names = {'Ip','rmag','zmag'};
C = [zeros(2,n); model.Cc(model.ind_zIp,:)/model.Ip0];

% operating point: equilibrium coil currents; the voltages are zero (descriptor without
% voltage entry): report the circuit voltages that would hold the equilibrium currents
Ia0 = model.Ia0(1:numel(coil_names));
current = struct();
for ii = 1:numel(coil_names)
    current.(coil_names{ii}) = Ia0(ii);
end
Vc0 = full(model.R(1:nc,1:nc))*model.Ic0(1:nc);
if max(abs(Vc0)) > 1e-6
    warning('write_vertical_circuit_model:voltage', ['Holding the equilibrium currents ', ...
        'needs circuit voltages up to %.3g V; the descriptor has zero operating-point ', ...
        'voltages, so these currents decay in the run.'], max(abs(Vc0)));
end
if any(model.Ic0(nc+1:end) ~= 0)
    warning('write_vertical_circuit_model:passive', ['The equilibrium has passive currents ', ...
        '(max |I| %.4g A): the plant reports only their deviations.'], max(abs(model.Ic0(nc+1:end))));
end

growth = max(real(eig(model.Ac)));
n_unstable = nnz(real(eig(model.Ac)) > 1e-9);
if isempty(shot)
    tok = regexp(char(equilibrium.uri),'(?<!\d)(\d{6})(?!\d)','tokens','once');
    shot = 0;
    if ~isempty(tok), shot = str2double(tok{1}); end
end
source_code = ['linear_models/utils/vertical_circuit (copy of plasma-tokamak-circuits, branch ', ...
    'vertical_instability, commit c02e88b)'];

d.description = sprintf(['vertical-instability model with circuit states around %s at %.10g s: ', ...
    'rigid vertical plasma displacement (Ip fixed, no radial motion), %d active circuits ', ...
    '(CS1U/CS1L in series, VS3U/VS3L in antiseries) and %d passive loops, the 14 power-supply ', ...
    'voltages as inputs; continuous time, dI/dt = A I + B U, deviations from the equilibrium ', ...
    'currents'], equilibrium.uri, equilibrium_time, nc, nl);
d.shot = shot;
d.time = equilibrium_time;
d.files = struct('A','A.npy','B','B.npy','C','C.npy');
d.states = [circuit_names, loop_names];
d.inputs = supply_names;
d.coil_states = coil_states;
d.passive_states = loop_names;
d.outputs = y_names;
d.coordinates = 'deviation';
d.operating_point = struct('coil_current', current);
d.units = struct('currents','A','voltages','V','rmag_zmag','m','Ip','A');
d.provenance = sprintf(['configure_vertical_instability_model of %s (A. Coutlis et al., Nucl. ', ...
    'Fusion 39 (1999); F. Pesamosca, PhD thesis, EPFL (2021), App. B); equilibrium %s at ', ...
    '%.10g s, passive currents: %s; written by write_vertical_circuit_model.m on %s'], ...
    source_code, equilibrium.uri, equilibrium_time, equilibrium.passive_source, ...
    char(datetime('today','Format','yyyy-MM-dd')));
d.notes = { ...
    ['Inputs: the 14 power-supply voltages, in the order the controller sends them; ', ...
     'Vc_CS1 = Vps_CS1U + Vps_CS1L, VS1 acts as -PF2 -PF3 +PF4 +PF5, Vc_VS3 = Vps_VS3.'], ...
    ['zmag = zIp/Ip0: rigid vertical displacement of the plasma current; the magnetic axis ', ...
     'moves by the same amount. The Ip and rmag rows are zero.'], ...
    sprintf(['%d unstable eigenvalue(s) of A, growth rate %.4g 1/s (e-folding time %.1f ms). ', ...
     'Operating point: the equilibrium coil currents, no voltages.'], n_unstable, growth, ...
     1e3/growth)};

if ~isfolder(out_dir)
    mkdir(out_dir);
end
write_npy(fullfile(out_dir,'A.npy'), model.Ac);
write_npy(fullfile(out_dir,'B.npy'), model.Bc);
write_npy(fullfile(out_dir,'C.npy'), C);
fid = fopen(fullfile(out_dir,'linear_model.json'),'w');
if fid < 0
    error('write_vertical_circuit_model:file','Cannot write %s.', ...
        fullfile(out_dir,'linear_model.json'));
end
fprintf(fid,'%s\n',jsonencode(d,'PrettyPrint',true));
fclose(fid);

% check: the files read back equal the model matrices
err = max([max(abs(read_npy_file(fullfile(out_dir,'A.npy')) - model.Ac),[],'all'), ...
    max(abs(read_npy_file(fullfile(out_dir,'B.npy')) - model.Bc),[],'all'), ...
    max(abs(read_npy_file(fullfile(out_dir,'C.npy')) - C),[],'all')]);
if err ~= 0
    error('write_vertical_circuit_model:check','The files read back differ from the model (%g).', err);
end
fprintf(['write_vertical_circuit_model: wrote %s\n  %d states (%d circuits + %d passive ', ...
    'loops), 14 power-supply inputs, outputs Ip, rmag, zmag\n  %d unstable eigenvalue(s), ', ...
    'growth rate %.4g 1/s; VS1 column = -PF2 -PF3 +PF4 +PF5 (checked)\n  use it with: ', ...
    'plant.linear_model_dir: "%s"\n'], out_dir, n, nc, nl, n_unstable, growth, out_dir);
end

function [value, opts] = take_option(opts, name, default)
% value of opts.(name) (removed from opts) or the default
value = default;
if isfield(opts, name)
    value = opts.(name);
    if ischar(default) || isstring(value)
        value = char(value);
    end
    opts = rmfield(opts, name);
end
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

function M = read_npy_file(filename)
% reads back what write_npy writes (float64, 2-D, Fortran order)
fid = fopen(filename,'r','ieee-le');
cleanup = onCleanup(@() fclose(fid));
fread(fid,8,'uint8');
header = char(fread(fid,fread(fid,1,'uint16'),'char')');
shape = sscanf(regexp(header,'\((\d+), (\d+)\)','match','once'),'(%d, %d)')';
M = reshape(fread(fid,prod(shape),'float64'),shape);
end
