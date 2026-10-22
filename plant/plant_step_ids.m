function [pf_active_out, pf_passive_out, x_new, y_out] = plant_step_ids(model, pf_active_in, x)
% one time step of a linear plant model (configure_linear_model) with IDS input and output
% model:        from configure_linear_model
% pf_active_in: one slice (homogeneous_time 1) with the 14 power-supply voltages u at t,
%               coil{k}.voltage = supply k (as the controller writes them, read by index)
% x:            state at t (model.n_states x 1)
% Step: du = u - op.voltage; x_new = A*x + B*du; the outputs at t+dt are
%   coil currents  op.coil_current + x_to_coil*x_new   (pf_active_out; the input
%                  voltages are echoed in the same slots, they are applied over [t, t+dt])
%   loop currents  x_to_loop*x_new                     (pf_passive_out)
%   observables    op.y + C_y*x_new                    (y_out: names, values, known)
% The states are deviations from the operating point op (zero for an absolute model).

if pf_active_in.ids_properties.homogeneous_time ~= 1 || numel(pf_active_in.time) ~= 1
    error('linear_model:slice', ['pf_active input must have homogeneous_time 1 and ', ...
        'exactly one time slice (homogeneous_time %d, %d times).'], ...
        pf_active_in.ids_properties.homogeneous_time, numel(pf_active_in.time));
end
t = pf_active_in.time;
if ~isequal(size(x),[model.n_states 1])
    error('linear_model:state','State must be a %dx1 vector.', model.n_states);
end

% the 14 power-supply voltages (model.input_names), by index, as the controller writes
% them: coil{k}.voltage carries power supply k
n_u = numel(model.input_names);
if numel(pf_active_in.coil) < n_u
    error('linear_model:voltage','pf_active has %d coils, %d power-supply voltages expected.', ...
        numel(pf_active_in.coil), n_u);
end
u = zeros(n_u,1);
for ii = 1:n_u
    v = pf_active_in.coil{ii}.voltage.data;
    if numel(v) ~= 1 || v == -9e40 || ~isfinite(v)
        error('linear_model:voltage','pf_active voltage %s: expected one value at t = %g, got %d.', ...
            model.input_names{ii}, t, numel(v));
    end
    u(ii) = v;
end

% state update (zero-order hold, model.A/B from c2d in configure_linear_model)
du = u - model.op.voltage;
x_new = model.A*x + model.B*du;
i_coil = model.op.coil_current + model.x_to_coil*x_new;
i_loop = model.x_to_loop*x_new;
y_out.names = model.y_names;
y_out.values = model.op.y + model.C_y*x_new;
y_out.known = model.op.y_known;

% outputs
label = ['linear model (' model.source.description ')'];
pf_active_out = ids_init('pf_active');
pf_active_out.ids_properties.homogeneous_time = 1;
pf_active_out.ids_properties.comment = [label ': coil currents'];
pf_active_out.time = t + model.dt;
pf_active_out.coil = ids_allocate('pf_active','coil',model.n_coils);
for ii = 1:model.n_coils
    pf_active_out.coil{ii}.name = model.coil_names{ii};
    if isfield(pf_active_out.coil{ii},'identifier')
        pf_active_out.coil{ii}.identifier = model.coil_names{ii};
    end
    pf_active_out.coil{ii}.current.data = i_coil(ii);
    pf_active_out.coil{ii}.voltage.data = u(ii);
end

pf_passive_out = ids_init('pf_passive');
pf_passive_out.ids_properties.homogeneous_time = 1;
pf_passive_out.ids_properties.comment = [label ': passive loop currents'];
pf_passive_out.time = t + model.dt;
pf_passive_out.loop = ids_allocate('pf_passive','loop',model.n_loops);
for ii = 1:model.n_loops
    pf_passive_out.loop{ii}.name = model.loop_names{ii};
    pf_passive_out.loop{ii}.current = i_loop(ii);
end

end
