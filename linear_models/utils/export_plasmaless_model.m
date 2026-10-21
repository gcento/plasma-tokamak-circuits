function export_plasmaless_model(plasmaless_repo, out_dir, md_uri)
% Export the plasmaless coil+vessel model (no plasma) of a machine description as the raw
% source linear_models/sources/plasmaless_model (linear_model.json + A, B, C .npy):
% 14 coils + passive loops, with the 14 coil voltages as inputs. The plant does not load
% it (its inputs are not the power supplies): linear_models/utils/
% transform_plasmaless_model.m reduces it to the ITER circuits and writes
% linear_models/plasmaless_model.
%   export_plasmaless_model(plasmaless_repo, out_dir)          % md: <repo>/data/md_dd4
%   export_plasmaless_model(plasmaless_repo, out_dir, md_uri)
% with out_dir = <repo>/linear_models/sources/plasmaless_model.
% Needs IMAS-MATLAB (reads the machine description) and the Control System Toolbox; run
% it on a compute node.
%
% The model is absolute (states = coil and passive currents), so the descriptor has
% coordinates 'absolute' and no operating point. After writing, A.npy and B.npy are read
% back and compared with the computed matrices.

plasmaless_repo = char(plasmaless_repo);
out_dir = char(out_dir);
if nargin < 3 || isempty(md_uri)
    md_uri = ['imas:hdf5?path=' fullfile(plasmaless_repo,'data','md_dd4')];
end
md_uri = char(md_uri);   % imas_open needs a char URI
addpath(fullfile(plasmaless_repo,'imas_model'));

% dt is not in the exported matrices: configure_plasmaless_model needs one (M and R do not
% depend on it)
dt = 0.002;
[em_coupling, pf_active, pf_passive] = load_plasmaless_md(md_uri);
native = configure_plasmaless_model(em_coupling, pf_active, pf_passive, dt);
n = native.n_states;
n_coils = native.n_coils;
% continuous matrices exactly as configure_plasmaless_model builds them
Ac = -native.M\native.R;
Bc = native.M\eye(n,n_coils);

if ~isfolder(out_dir)
    mkdir(out_dir);
end
write_npy(fullfile(out_dir,'A.npy'), Ac);
write_npy(fullfile(out_dir,'B.npy'), Bc);
write_npy(fullfile(out_dir,'C.npy'), zeros(0,n));

d.description = ['plasmaless coil+vessel circuit model (no plasma), ', num2str(n_coils), ...
    ' coils + ', num2str(native.n_loops), ' passive loops: dI/dt = A I + B U, ', ...
    'absolute currents and voltages'];
d.shot = 0;
d.time = 0;
d.files = struct('A','A.npy','B','B.npy','C','C.npy');
d.states = native.state_names;
d.inputs = native.coil_names;
d.coil_states = struct();
for ii = 1:n_coils
    d.coil_states.(native.coil_names{ii}) = native.coil_names{ii};
end
d.passive_states = native.loop_names;
d.outputs = {};
d.coordinates = 'absolute';
d.units = struct('currents','A','voltages','V');
d.provenance = sprintf(['machine description %s; A = -M\\R, B = M\\eye(n,n_coils) as ', ...
    'configure_plasmaless_model.m; written by linear_models/utils/export_plasmaless_model.m on %s. ', ...
    'Raw source: the inputs are the 14 coil voltages, so the plant does not load it; ', ...
    'linear_models/utils/transform_plasmaless_model.m writes its circuit form ', ...
    'linear_models/plasmaless_model'], ...
    md_uri, char(datetime('today','Format','yyyy-MM-dd')));
d.notes = {['Not tied to a pulse: shot and time are 0. Absolute model (no operating point): ', ...
    'the actor starts it from the F_INIT coil currents at t0 (passive currents 0) and ', ...
    'ignores the setting operating_point.']};
fid = fopen(fullfile(out_dir,'linear_model.json'),'w');
fprintf(fid,'%s\n',jsonencode(d,'PrettyPrint',true));
fclose(fid);

% check: the matrix files read back as transform_plasmaless_model.m reads them
if ~isequal(read_npy(fullfile(out_dir,'A.npy')), Ac) || ~isequal(read_npy(fullfile(out_dir,'B.npy')), Bc)
    error('export_plasmaless_model:check','A.npy or B.npy differ from the computed matrices.');
end
fprintf(['export_plasmaless_model: wrote %s (%d states, %d coil-voltage inputs); next: ', ...
    'linear_models/utils/transform_plasmaless_model.m\n'], out_dir, n, n_coils);
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
fwrite(fid,M(:),'float64');
fclose(fid);
end
