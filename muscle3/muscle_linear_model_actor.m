function muscle_linear_model_actor()
% MUSCLE3 actor that runs one linear plant model of linear_models/ in place of NICE
% direct evolutive (nice_evo_rd), for controller tests in PDS (workflow
% linear_plant_controller). Run with
%   matlab -batch "pyenv('Version',<python with muscle3 0.10>,'ExecutionMode','InProcess'); ...
%                  addpath(<repo>/muscle3); muscle_linear_model_actor"
% Ports:
%   F_INIT equilibrium_f_init  reference trajectory; its timestamp is the start time t0
%          pf_active_f_init    coil currents of the scenario
%   S      pf_active_s         the 14 power-supply voltages, coil{k}.voltage = supply k
%                              (CS3U CS2U CS1U CS1L CS2L CS3L PF1..PF5 VS1 PF6 VS3)
%   O_I    equilibrium_o_i     reference ip and geometric axis + the model outputs
%          pf_active_o_i       coil currents
%          pf_passive_o_i      passive currents
% Settings: dt, t_interval (= dt), t_end, linear_model_dir (mandatory); md_uri (default
% data/md_dd4 of this repository), operating_point ('descriptor' or 'f_init').
% Operating point, initial state, timing (as NICE) and assumptions: muscle3/README.md.
% Log: one line for the model, one for the operating point, one at the end.

% plant/ next to this folder: model loader and step function (not a setting)
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo_root,'plant'));

ports = py.dict();
ports{py.getattr(@py.ymmsl.Operator,"F_INIT")} = py.list({"equilibrium_f_init", "pf_active_f_init"});
ports{py.getattr(@py.ymmsl.Operator,"S")} = py.list({"pf_active_s"});
ports{py.getattr(@py.ymmsl.Operator,"O_I")} = py.list({"equilibrium_o_i", "pf_active_o_i", "pf_passive_o_i"});
flags = py.getattr(@py.libmuscle.InstanceFlags,"KEEPS_NO_STATE_FOR_NEXT_USE");
instance = py.libmuscle.Instance(ports, flags);

try
    while instance.reuse_instance()
        run_one(instance);
    end
catch ME
    fprintf(2,'muscle_linear_model_actor: %s\n', getReport(ME,'extended','hyperlinks','off'));
    instance.error_shutdown(ME.message);
    rethrow(ME);
end

end

function run_one(instance)
% one reuse of the instance: F_INIT, model build, time loop

% settings
dt = double(instance.get_setting('dt'));
t_interval = double(instance.get_setting('t_interval'));
t_end = double(instance.get_setting('t_end'));
% default machine description: data/md_dd4 shipped in this repository
md_uri = char(instance.get_setting('md_uri', pyargs('default', ['imas:hdf5?path=' ...
    fullfile(fileparts(fileparts(mfilename('fullpath'))),'data','md_dd4')])));
% folder of the linear model (mandatory)
linear_model_dir = char(instance.get_setting('linear_model_dir', pyargs('default', '')));
if isempty(linear_model_dir)
    error('linear_model_actor:settings', ['linear_model_dir is empty: set it to a model ', ...
        'folder of linear_models/ (linear_models/plasmaless_model for the plasma-free ', ...
        'circuit model).']);
end
% which source wins for the operating point of a linear model
op_source = char(instance.get_setting('operating_point', pyargs('default', 'descriptor')));
if ~any(strcmp(op_source, {'descriptor','f_init'}))
    error('linear_model_actor:settings','operating_point must be ''descriptor'' or ''f_init'' (got ''%s'').', op_source);
end
if ~(dt > 0)
    error('linear_model_actor:settings','dt must be > 0 (got %g).', dt);
end
if t_interval ~= dt
    % the model output is only defined at t0 + k*dt; an exchange period other
    % than dt would need sub-stepping that NICE does but this actor does not
    error('linear_model_actor:settings','t_interval (%g) must equal dt (%g).', t_interval, dt);
end

% F_INIT
msg_eq = instance.receive("equilibrium_f_init");
t0 = double(msg_eq.timestamp);
eq_ref = imas_deserialize(uint8(msg_eq.data),'equilibrium');
msg_pfa = instance.receive("pf_active_f_init");
pfa_init = imas_deserialize(uint8(msg_pfa.data),'pf_active');

[eq_time, ip_ref, r_ref, z_ref] = reference_trajectory(eq_ref);
if ~(t_end > t0)
    error('linear_model_actor:t_end','t_end (%.10g) must be > t0 (%.10g).', t_end, t0);
end
if t0 < eq_time(1) || t_end > eq_time(end)
    error('linear_model_actor:t_end', ['[t0, t_end] = [%.10g, %.10g] is not inside the ', ...
        'F_INIT equilibrium time range [%.10g, %.10g].'], t0, t_end, eq_time(1), eq_time(end));
end

% machine description: only pf_active, for the coil names and order
idx = imas_open(md_uri,40);
pf_active_md = ids_get(idx,'pf_active');
imas_close(idx);

% The controller reads and writes the coils by index, in the order of the F_INIT
% pf_active; the model works in machine-description order. Any difference between the
% two would send currents to the wrong coils: error.
names_init = cellfun(@short_name, pfa_init.coil(:)', 'UniformOutput', false);
if ~isequal(names_init, model_coil_names(pf_active_md))
    error('linear_model_actor:coilOrder', ['F_INIT pf_active coils {%s} differ from the ', ...
        'machine-description coils {%s} (names or order).'], ...
        strjoin(names_init,','), strjoin(model_coil_names(pf_active_md),','));
end
coil_names = model_coil_names(pf_active_md);

pfa_t0 = closest_slice(pfa_init, t0, coil_names);
% x = 0 (deviations from the operating point), or the initial state of an absolute model
[model, x] = configure_linear_model(linear_model_dir, dt, pf_active_md);
fprintf('muscle_linear_model_actor: %s, %d states, t0 = %.10g, t_end = %.10g, dt = %g\n', ...
    linear_model_dir, model.n_states, t0, t_end, dt);
if strcmp(model.coordinates,'absolute')
    % no operating point: start from the F_INIT coil currents at t0 projected on the
    % circuit states (error if they violate the circuits), passive states 0.
    % ASSUMPTION: the scenario data has no passive-structure currents.
    i0 = zeros(model.n_coils,1);
    for ii = 1:model.n_coils
        i0(ii) = check_value(pfa_t0.coil{ii}.current.data, ...
            ['pf_active coil ' model.coil_names{ii} ' current'], t0);
    end
    x = linear_initial_state(model, i0);
    fprintf('muscle_linear_model_actor: absolute model, initial state from the F_INIT coil currents\n');
else
    model = linear_operating_point(model, pfa_t0, eq_ref, eq_time, t0, op_source);
end
eq_comment = sprintf(['linear model (%s): %s from the model, ', ...
    'boundary geometric axis = reference pass-through'], model.source.description, ...
    strjoin(model.y_names(model.op.y_known),'/'));

% time loop
pfa_in = pfa_t0;
% first step: the operating-point power-supply voltages (du = 0); the F_INIT coil
% voltages are scenario coil voltages, not power-supply voltages
for ii = 1:numel(model.op.voltage)
    pfa_in.coil{ii}.voltage.data = model.op.voltage(ii);
end
t = t0;
n_sent = 0;
while true
    [pfa_out, pfp_out, x, y_out] = plant_step_ids(model, pfa_in, x);
    t = pfa_out.time;      % = t + dt
    eq_out = reference_equilibrium(t, eq_time, ip_ref, r_ref, z_ref);
    % known observables of the model (absolute values): a known Ip overrides the
    % pass-through ip; magnetic axis, X-point and gaps are added;
    % boundary.geometric_axis stays the reference (the model does not give it)
    eq_out = plant_observables_to_equilibrium(eq_out, y_out);
    eq_out.ids_properties.comment = eq_comment;

    more_after_this = t + t_interval <= t_end + 1e-9;   % as cc:184
    if more_after_this
        t_next = py.float(t + t_interval);
    else
        t_next = py.None;
    end
    send_ids(instance, "equilibrium_o_i", eq_out, 'equilibrium', t, t_next);
    send_ids(instance, "pf_active_o_i", pfa_out, 'pf_active', t, t_next);
    send_ids(instance, "pf_passive_o_i", pfp_out, 'pf_passive', t, t_next);
    n_sent = n_sent + 1;
    if ~more_after_this
        break
    end

    % next input; the stop decision comes from t_end only (the controller always
    % sends next_timestamp None)
    msg = instance.receive("pf_active_s");
    pfa_in = imas_deserialize(uint8(msg.data),'pf_active');
    t_msg = double(msg.timestamp);
    % the controller stamps its command with the time of the measurement it answers;
    % 1e-6 s: well below dt, above the round-off of the accumulated time
    if abs(t_msg - t) > 1e-6 || numel(pfa_in.time) ~= 1 || abs(pfa_in.time - t) > 1e-6
        error('linear_model_actor:time', ['pf_active_s at timestamp %.10g (IDS time %s) does not ', ...
            'answer the measurement sent at t = %.10g.'], t_msg, mat2str(pfa_in.time), t);
    end
    pfa_in.time = t;       % keep this actor's own clock (differences < 1e-6 s)
    % the controller sends no coil names: set them by index (F_INIT order, checked
    % above to be the machine-description order), as NICE does
    if numel(pfa_in.coil) ~= numel(coil_names)
        error('linear_model_actor:coils','pf_active_s has %d coils, the model %d.', ...
            numel(pfa_in.coil), numel(coil_names));
    end
    for ii = 1:numel(pfa_in.coil)
        pfa_in.coil{ii}.name = coil_names{ii};
    end
end
fprintf('muscle_linear_model_actor: %d exchanges, last t = %.10g\n', n_sent, t);

end

function [eq_time, ip_ref, r_ref, z_ref] = reference_trajectory(eq)
% ip and boundary geometric axis of every F_INIT equilibrium time slice
if eq.ids_properties.homogeneous_time ~= 1
    error('linear_model_actor:equilibrium','F_INIT equilibrium must have homogeneous_time 1 (got %d).', ...
        eq.ids_properties.homogeneous_time);
end
eq_time = eq.time(:);
n = numel(eq.time_slice);
if n ~= numel(eq_time) || n == 0
    error('linear_model_actor:equilibrium','F_INIT equilibrium has %d time slices and %d times.', ...
        n, numel(eq_time));
end
if any(diff(eq_time) <= 0)
    error('linear_model_actor:equilibrium','F_INIT equilibrium times are not increasing.');
end
ip_ref = zeros(n,1); r_ref = zeros(n,1); z_ref = zeros(n,1);
for ii = 1:n
    s = eq.time_slice{ii};
    ip_ref(ii) = check_value(s.global_quantities.ip, 'equilibrium global_quantities.ip', eq_time(ii));
    r_ref(ii) = check_value(s.boundary.geometric_axis.r, 'equilibrium boundary.geometric_axis.r', eq_time(ii));
    z_ref(ii) = check_value(s.boundary.geometric_axis.z, 'equilibrium boundary.geometric_axis.z', eq_time(ii));
end
end

function v = check_value(v, label, t)
if numel(v) ~= 1 || v == -9e40 || ~isfinite(v)
    error('linear_model_actor:f_init','F_INIT %s missing at t = %.10g.', label, t);
end
end

function eq = reference_equilibrium(t, eq_time, ip_ref, r_ref, z_ref)
% one-slice equilibrium: ip and geometric axis of the reference trajectory,
% linearly interpolated at t (reference pass-through); a known Ip of the model replaces
% ip afterwards (plant_observables_to_equilibrium).
if numel(eq_time) == 1
    if t ~= eq_time
        error('linear_model_actor:equilibrium','t = %.10g outside the one-slice F_INIT equilibrium.', t);
    end
    vals = [ip_ref r_ref z_ref];
else
    vals = interp1(eq_time, [ip_ref r_ref z_ref], t, 'linear');
end
if any(~isfinite(vals))
    error('linear_model_actor:equilibrium','t = %.10g outside the F_INIT equilibrium time range.', t);
end
eq = ids_init('equilibrium');
eq.ids_properties.homogeneous_time = 1;
eq.ids_properties.comment = 'reference pass-through of ip and geometric axis';
eq.time = t;
eq.time_slice = ids_allocate('equilibrium','time_slice',1);
eq.time_slice{1}.time = t;
eq.time_slice{1}.global_quantities.ip = vals(1);
eq.time_slice{1}.boundary.geometric_axis.r = vals(2);
eq.time_slice{1}.boundary.geometric_axis.z = vals(3);
eq.code.output_flag = 0;
end

function model = linear_operating_point(model, pfa_t0, eq_ref, eq_time, t0, op_source)
% operating point of a linear model: descriptor entries (model.op, *_given /
% y_known flags) and the F_INIT data at t0 (coil currents of the pf_active slice
% pfa_t0 in model coil order, observables of the equilibrium slice closest to t0).
% op_source 'descriptor': descriptor entries win, F_INIT fills the rest; 'f_init':
% F_INIT values win, descriptor entries fill the rest. The power-supply voltages
% come from the descriptor only (zero when absent): the F_INIT coil voltages are
% scenario coil voltages, not power-supply voltages.
f_init_wins = strcmp(op_source,'f_init');
n_cur = 0;    % coil currents taken from F_INIT
for ii = 1:model.n_coils
    c = pfa_t0.coil{ii};
    if f_init_wins || ~model.op.coil_current_given(ii)
        model.op.coil_current(ii) = check_value(c.current.data, ...
            ['pf_active coil ' model.coil_names{ii} ' current'], t0);
        n_cur = n_cur + 1;
    end
end
[~, k] = min(abs(eq_time - t0));
[vals, known] = plant_observables_from_equilibrium(eq_ref.time_slice{k}, model.y_names);
if f_init_wins
    from_eq = known;
else
    from_eq = known & ~model.op.y_known;
end
model.op.y(from_eq) = vals(from_eq);
model.op.y_known = model.op.y_known | known;
fprintf(['muscle_linear_model_actor: operating_point = %s: %d of %d coil currents from ', ...
    'F_INIT; outputs sent {%s}, not sent {%s}\n'], op_source, n_cur, model.n_coils, ...
    strjoin(model.y_names(model.op.y_known),','), strjoin(model.y_names(~model.op.y_known),','));
end

function names = model_coil_names(pf_active_md)
% short coil names of the machine description, in its order (the model order)
names = cellfun(@short_name, pf_active_md.coil(:)', 'UniformOutput', false);
end

function pfa = closest_slice(pfa_all, t, coil_names)
% one-slice pf_active at the sample of pfa_all closest to t (as NICE's
% getSlice(..., CLOSEST_SAMPLE), cc:156; equidistant samples: the earlier one)
% with coil currents and voltages, coils named coil_names (by index)
if pfa_all.ids_properties.homogeneous_time ~= 1
    error('linear_model_actor:pf_active','F_INIT pf_active must have homogeneous_time 1 (got %d).', ...
        pfa_all.ids_properties.homogeneous_time);
end
tt = pfa_all.time(:);
if isempty(tt)
    error('linear_model_actor:pf_active','F_INIT pf_active has no time.');
end
[~, k] = min(abs(tt - t));
pfa = ids_init('pf_active');
pfa.ids_properties.homogeneous_time = 1;
pfa.time = t;
n = numel(pfa_all.coil);
pfa.coil = ids_allocate('pf_active','coil',n);
for ii = 1:n
    c = pfa_all.coil{ii};
    if numel(c.current.data) ~= numel(tt) || numel(c.voltage.data) ~= numel(tt)
        error('linear_model_actor:pf_active', ['F_INIT pf_active coil %d: %d currents and %d ', ...
            'voltages for %d times.'], ii, numel(c.current.data), numel(c.voltage.data), numel(tt));
    end
    pfa.coil{ii}.name = coil_names{ii};
    pfa.coil{ii}.current.data = c.current.data(k);
    pfa.coil{ii}.voltage.data = c.voltage.data(k);
end
if abs(tt(k) - t) > 1e-9
    fprintf('muscle_linear_model_actor: no F_INIT pf_active sample at t0, using t = %.10g\n', tt(k));
end
end

function send_ids(instance, port, ids, ids_name, t, t_next)
% serialization as KCURR_RZIp/muscle_NICE_input.m:53-57
serialized = imas_serialize(ids, ids_name);
s = sprintf('%d ' ,uint8(serialized));
x = py.numpy.fromstring(s, py.numpy.int8, int8(-1), char(' '));
serialized = py.bytes(x);
msg = py.libmuscle.Message(t, t_next, serialized);
instance.send(port, msg);
end

function key = short_name(coil)
% short coil key, the rule shared with plant/ (plant_coil_key.m)
key = plant_coil_key(coil);
end
