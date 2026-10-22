function [model, x0] = configure_linear_model(model_dir, dt, pf_active_md)
% Load a model folder of linear_models/ (linear_model.json + A, B, C) and build the
% discrete plant model of  dx/dt = A x + B u,  y = C x  (zero-order hold with step dt):
% x the circuit and passive states, u the 14 power-supply voltages, y the plasma outputs.
% Descriptor fields, circuit contract, sign convention and errors: linear_models/README.md.
%   [model, x0] = configure_linear_model(model_dir, dt, pf_active_md)
% pf_active_md: machine-description pf_active; its coils (names and order) are
% model.coil_names.
% x0 = 0: with coordinates 'deviation' (default) x, u and y are deviations from the
% operating point model.op; with 'absolute' there is no operating point and the caller
% starts instead from the coil currents (linear_initial_state).
% model: kind, coordinates, sysd, A, B, C, D (discrete), dt, n_states, state_names,
% coil_names, n_coils, x_to_coil (coil currents = x_to_coil*x), loop_names, n_loops,
% x_to_loop, y_names, C_y, input_names, op, source, convention; A_c, B_c, C_c
% (continuous, IMAS convention); A_file, B_file, C_file (as read from the files).

if nargin ~= 3
    error('linear_model:nargin','Use configure_linear_model(model_dir, dt, pf_active_md).');
end
desc_file = fullfile(model_dir,'linear_model.json');
if ~isfile(desc_file)
    error('linear_model:descriptor','Descriptor %s not found.', desc_file);
end
d = jsondecode(fileread(desc_file));
required = {'description','shot','time','files','states','inputs','coil_states', ...
    'passive_states','outputs'};
for ii = 1:numel(required)
    if ~isfield(d,required{ii})
        error('linear_model:descriptor','%s: field ''%s'' missing.', desc_file, required{ii});
    end
end
for f = {'A','B','C'}
    if ~isfield(d.files,f{1})
        error('linear_model:descriptor','%s: files.%s missing.', desc_file, f{1});
    end
end
states = name_list(d.states,'states',desc_file);
inputs = name_list(d.inputs,'inputs',desc_file);
passive = name_list(d.passive_states,'passive_states',desc_file);
outputs = name_list(d.outputs,'outputs',desc_file);
n = numel(states); m = numel(inputs); p = numel(outputs);
if ~isstruct(d.coil_states)
    error('linear_model:descriptor','%s: coil_states must be an object coil -> state.', desc_file);
end

% outputs vocabulary (before the convention block, which refers to the outputs)
for ii = 1:p
    if isempty(regexp(outputs{ii},'^(rmag|zmag|rX|zX|Ip|Dist2ndSep|gap[1-9][0-9]*)$','once'))
        error('linear_model:observable', ['%s: output ''%s'' is not in the vocabulary ', ...
            '{rmag, zmag, rX, zX, Ip, gap<k>, Dist2ndSep}.'], desc_file, outputs{ii});
    end
end

% matrices
A = load_matrix(fullfile(model_dir,d.files.A));
B = load_matrix(fullfile(model_dir,d.files.B));
C = load_matrix(fullfile(model_dir,d.files.C));
if ~isequal(size(A),[n n]) || ~isequal(size(B),[n m]) || ~isequal(size(C),[p n])
    error('linear_model:size', ['%s: A %s, B %s, C %s do not match %d states, %d inputs, ', ...
        '%d outputs.'], desc_file, mat2str(size(A)), mat2str(size(B)), mat2str(size(C)), n, m, p);
end
if ~all(isfinite(A(:))) || ~all(isfinite(B(:))) || ~all(isfinite(C(:)))
    error('linear_model:size','%s: the matrices contain non-finite values.', desc_file);
end

% sign convention: from the frame of the files to IMAS (linear_models/README.md)
conv = sign_convention(d, outputs, desc_file);
A_file = A; B_file = B; C_file = C;
B = (conv.input_sign/conv.state_sign)*B;
C = (conv.state_sign./conv.output_sign).*C;   % row k scaled by state_sign/output_sign(k)

% coils: machine-description names and order
n_coils = numel(pf_active_md.coil);
coil_names = cell(1,n_coils);
for ii = 1:n_coils
    coil_names{ii} = plant_coil_key(pf_active_md.coil{ii});
end
if numel(unique(coil_names)) ~= n_coils
    error('linear_model:coilName','Duplicate pf_active coil names.');
end
coil_keys = fieldnames(d.coil_states)';
if ~isequal(sort(coil_keys),sort(coil_names))
    error('linear_model:descriptor','%s: coil_states keys {%s} are not the machine-description coils {%s}.', ...
        desc_file, strjoin(coil_keys,','), strjoin(coil_names,','));
end

% inputs: the power-supply voltages in the order the controller sends them
supply_names = {'CS3U','CS2U','CS1U','CS1L','CS2L','CS3L','PF1','PF2','PF3','PF4', ...
    'PF5','VS1','PF6','VS3'};
if ~isequal(inputs, supply_names)
    error('linear_model:descriptor', ['%s: inputs {%s} must be the power supplies {%s}, in ', ...
        'this order (the order the controller sends them).'], desc_file, ...
        strjoin(inputs,','), strjoin(supply_names,','));
end
if isfield(d,'input_map')
    error('linear_model:descriptor', ['%s: input_map is no longer used: the inputs are the ', ...
        'power-supply voltages, applied to B as they are.'], desc_file);
end

% states: every state is a coil state or a passive state, not both
x_to_coil = zeros(n_coils,n);
is_coil_state = false(1,n);
for ii = 1:n_coils
    s = d.coil_states.(coil_names{ii});
    if ~ischar(s) && ~isstring(s)
        error('linear_model:descriptor','%s: coil_states.%s must be a state name.', desc_file, coil_names{ii});
    end
    s = char(s);
    k = find(strcmp(states,s),1);
    sgn = 1;
    if isempty(k) && startsWith(s,'-')
        % '-name': minus the state name (coil in antiseries, e.g. VS3U = -VS3); a
        % state literally called '-name' (as in the NICE models) takes precedence
        k = find(strcmp(states,s(2:end)),1);
        sgn = -1;
    end
    if isempty(k)
        error('linear_model:descriptor','%s: coil_states.%s = ''%s'' is not a state.', ...
            desc_file, coil_names{ii}, s);
    end
    x_to_coil(ii,k) = sgn;
    is_coil_state(k) = true;
end
is_passive = ismember(states,passive);
missing = setdiff(passive,states);
if ~isempty(missing)
    error('linear_model:descriptor','%s: passive_states {%s} are not states.', desc_file, strjoin(missing,','));
end
both = states(is_coil_state & is_passive);
if ~isempty(both)
    error('linear_model:descriptor','%s: states {%s} are both coil and passive states.', desc_file, strjoin(both,','));
end
none = states(~is_coil_state & ~is_passive);
if ~isempty(none)
    error('linear_model:descriptor','%s: states {%s} are neither coil nor passive states.', desc_file, strjoin(none,','));
end
% ITER circuits, as in NICE: CS1U = CS1L (series), VS3U = -VS3L (antiseries), every
% other coil on a state of its own
pair_coils = {'CS1U','CS1L','VS3U','VS3L'};
row = @(c) x_to_coil(strcmp(coil_names,c),:);
n_share = sum(x_to_coil ~= 0,1);                   % coils on each state
own_ok = all(cellfun(@(c) n_share(row(c) ~= 0) == 1, setdiff(coil_names,pair_coils,'stable')));
if ~all(ismember(pair_coils,coil_names)) || ~own_ok || ...
        ~isequal(row('CS1U'),row('CS1L')) || ~isequal(row('VS3U'),-row('VS3L'))
    map = cellfun(@(c) sprintf('%s -> %s', c, char(d.coil_states.(c))), ...
        intersect(pair_coils,coil_names,'stable'), 'UniformOutput', false);
    error('linear_model:circuits', ['%s: coil_states do not follow the ITER circuits: ', ...
        'CS1U and CS1L on one state with the same sign (series), VS3U and VS3L on one ', ...
        'state with opposite signs (antiseries, e.g. "VS3" and "-VS3"), every other coil ', ...
        'on a state of its own. Got %s; states shared by several coils: {%s}.'], ...
        desc_file, strjoin(map,', '), strjoin(states(n_share > 1),','));
end
% coordinates: deviations from the operating point (default) or absolute values
coordinates = 'deviation';
if isfield(d,'coordinates')
    coordinates = char(d.coordinates);
    if ~any(strcmp(coordinates,{'deviation','absolute'}))
        error('linear_model:descriptor','%s: coordinates must be ''deviation'' or ''absolute'' (got ''%s'').', ...
            desc_file, coordinates);
    end
end
if strcmp(coordinates,'absolute')
    if isfield(d,'operating_point')
        error('linear_model:descriptor','%s: an absolute model has no operating_point.', desc_file);
    end
end

loop_names = states(is_passive);
n_loops = numel(loop_names);
x_to_loop = zeros(n_loops,n);
x_to_loop(:,is_passive) = eye(n_loops);

% discrete model (zero-order hold)
ev = eig(A);
if any(real(ev) > 1e-9)
    warning('linear_model:unstable','%s: %d eigenvalue(s) of A with real part > 1e-9 (max %.3g).', ...
        desc_file, nnz(real(ev) > 1e-9), max(real(ev)));
end
sys = ss(A,B,C,zeros(p,m));
sysd = c2d(sys,dt);
if any(sysd.D(:) ~= 0)
    error('linear_model:size','%s: the discrete D matrix is not zero.', desc_file);
end

model.kind = 'linear';
model.coordinates = coordinates;
model.sysd = sysd;
model.A = sysd.A;
model.B = sysd.B;
model.C = sysd.C;
model.D = sysd.D;
model.dt = dt;
model.n_states = n;
model.state_names = states;
model.coil_names = coil_names;
model.n_coils = n_coils;
model.x_to_coil = x_to_coil;
model.loop_names = loop_names;
model.n_loops = n_loops;
model.x_to_loop = x_to_loop;
model.y_names = outputs;
model.C_y = sysd.C;
model.input_names = inputs;
model.A_c = A;
model.B_c = B;
model.C_c = C;
model.convention = conv;
model.A_file = A_file;
model.B_file = B_file;
model.C_file = C_file;
model.source = struct('model_dir',model_dir,'files',d.files,'shot',d.shot,'time',d.time, ...
    'description',d.description);

% operating point: zeros, overridden by the descriptor (currents per coil,
% voltages per power supply)
op.coil_current = zeros(n_coils,1);
op.coil_current_given = false(n_coils,1);
op.voltage = zeros(m,1);
op.voltage_given = false(m,1);
op.y = zeros(p,1);
op.y_known = false(p,1);
if isfield(d,'operating_point') && isstruct(d.operating_point)
    o = d.operating_point;
    if isfield(o,'coil_current')
        [op.coil_current, op.coil_current_given] = op_values(o.coil_current,coil_names, ...
            op.coil_current,'operating_point.coil_current',desc_file);
    end
    if isfield(o,'voltage')
        [op.voltage, op.voltage_given] = op_values(o.voltage,inputs, ...
            op.voltage,'operating_point.voltage',desc_file);
    end
    if isfield(o,'outputs')
        [op.y, op.y_known] = op_values(o.outputs,outputs,op.y,'operating_point.outputs',desc_file);
    end
end
model.op = op;

x0 = zeros(n,1);

end

function names = name_list(v, label, desc_file)
% cellstr row of unique names from a jsondecode'd JSON array of strings
if isempty(v)
    names = {};
elseif ischar(v) || isstring(v)
    names = {char(v)};
elseif iscellstr(v)
    names = v(:)';
else
    error('linear_model:descriptor','%s: %s must be a list of names.', desc_file, label);
end
if numel(unique(names)) ~= numel(names)
    error('linear_model:descriptor','%s: %s has duplicate names.', desc_file, label);
end
end

function conv = sign_convention(d, outputs, desc_file)
% convention block: state_sign, input_sign (+1 or -1), output_sign (output ->
% +1 or -1); everything absent defaults to +1 (files in the wrapper's frame)
p = numel(outputs);
conv.state_sign = 1;
conv.input_sign = 1;
conv.output_sign = ones(p,1);
if ~isfield(d,'convention')
    return
end
c = d.convention;
if ~isstruct(c)
    error('linear_model:descriptor','%s: convention must be an object.', desc_file);
end
unknown = setdiff(fieldnames(c)', {'state_sign','input_sign','output_sign'});
if ~isempty(unknown)
    error('linear_model:descriptor','%s: convention has unknown field(s) {%s}.', ...
        desc_file, strjoin(unknown,','));
end
if isfield(c,'state_sign')
    conv.state_sign = sign_value(c.state_sign,'convention.state_sign',desc_file);
end
if isfield(c,'input_sign')
    conv.input_sign = sign_value(c.input_sign,'convention.input_sign',desc_file);
end
if isfield(c,'output_sign')
    if ~isstruct(c.output_sign)
        error('linear_model:descriptor','%s: convention.output_sign must be an object output -> +1/-1.', ...
            desc_file);
    end
    names = fieldnames(c.output_sign)';
    for ii = 1:numel(names)
        k = find(strcmp(outputs,names{ii}),1);
        if isempty(k)
            error('linear_model:descriptor','%s: convention.output_sign.%s is not an output.', ...
                desc_file, names{ii});
        end
        conv.output_sign(k) = sign_value(c.output_sign.(names{ii}), ...
            ['convention.output_sign.' names{ii}], desc_file);
    end
end
end

function s = sign_value(v, label, desc_file)
% a convention entry: the number +1 or -1
if ~isnumeric(v) || ~isscalar(v) || ~(v == 1 || v == -1)
    error('linear_model:descriptor','%s: %s must be +1 or -1.', desc_file, label);
end
s = double(v);
end

function M = load_matrix(filename)
if ~isfile(filename)
    error('linear_model:file','Matrix file %s not found.', filename);
end
[~, ~, ext] = fileparts(filename);
switch lower(ext)
    case '.npy'
        M = read_npy(filename);
    case '.txt'
        M = readmatrix(filename,'FileType','text');
    otherwise
        error('linear_model:descriptor','%s: matrix files must be .npy or .txt.', filename);
end
end

function [vals, given] = op_values(o, names, vals, label, desc_file)
% apply the key -> value entries of a descriptor object to vals (indexed by names)
if ~isstruct(o)
    error('linear_model:descriptor','%s: %s must be an object name -> value.', desc_file, label);
end
given = false(size(vals));
keys = fieldnames(o)';
for ii = 1:numel(keys)
    k = find(strcmp(names,keys{ii}),1);
    v = o.(keys{ii});
    if isempty(k)
        error('linear_model:descriptor','%s: %s.%s is not a known name.', desc_file, label, keys{ii});
    end
    if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v)
        error('linear_model:descriptor','%s: %s.%s must be a finite number.', desc_file, label, keys{ii});
    end
    vals(k) = v;
    given(k) = true;
end
end
