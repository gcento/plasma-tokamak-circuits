function res = supply_sequence_response(model_dir, seq, pf_active_md)
% Open-loop response of one model folder (linear_models/<name>) to voltages on the power
% supplies, one supply after the other: supply i of seq.supplies starts at
% t_on + (i-1)*spacing, every other supply stays at 0 V. With one supply it is the
% single-channel test. Computed as the plant computes it (configure_linear_model, then
% x(k+1) = A x(k) + B u(k) from x = 0): every result is a deviation from the operating
% point of the model. Launcher: utils/run_supply_sequence.m.
%   res = supply_sequence_response(model_dir, seq)
%   res = supply_sequence_response(model_dir, seq, pf_active_md)
% seq fields [default]:
%   supplies   power supplies in pulse order [all 14]
%   waveform   'bipolar' (+A for width/2, then -A for width/2), 'pulse' (+A for width)
%              or 'step' (+A until the end) ['bipolar']
%   amplitude  [V] one value or one per supply [10]
%   t_on       [s] start of the first supply [0.01]
%   spacing    [s] between the starts of two supplies [2]
%   width      [s] duration of pulse and bipolar [0.1]
%   t_end      [s] end of the run [t_on + n*spacing]
%   dt         [s] time step, the plant's dt [1e-3]
%   i_limit    [A] stop when a coil current exceeds it (unstable models) [1e6]
% The voltages are defined on the sample index k (t_k = k*dt, u_k applied over
% [t_k, t_k + dt]) with half-up rounding, exactly as muscle3/voltage_driver.py.
% res: label, model_dir, seq, time, supply_names, u (nt x 14), pulse_start, coil_names,
%   Ia (nt x 14 coil currents), circuit_names, Ic (nt x 12 circuit currents), y_names,
%   Y (outputs), n_states, growth (largest real part of the eigenvalues of A), stopped_at.
% Circuit currents as the controller measures them (utils/circuit_currents.m).

if nargin < 3 || isempty(pf_active_md)
    [~, pf_active_md] = load_machine_description();
end
seq = check_sequence(seq);
model_dir = char(model_dir);
w = warning('off','linear_model:unstable');        % vertical models: unstable by design
restore = onCleanup(@() warning(w));
model = configure_linear_model(model_dir, seq.dt, pf_active_md);
[known, j] = ismember(seq.supplies, model.input_names);
if ~all(known)
    error('supply_sequence_response:input', 'Unknown power supplies {%s}: use %s.', ...
        strjoin(seq.supplies(~known),','), strjoin(model.input_names,' '));
end

% voltages on the sample index
nt = round(seq.t_end/seq.dt) + 1;
k = (0:nt-1)';
k_start = floor((seq.t_on + (0:numel(j)-1)'*seq.spacing)/seq.dt + 0.5);
n_w = floor(seq.width/seq.dt + 0.5);
n_h = floor(n_w/2 + 0.5);
u = zeros(nt, numel(model.input_names));
for i = 1:numel(j)
    a = seq.amplitude(i);
    on = k >= k_start(i);
    switch seq.waveform
        case 'step'
            ui = a*on;
        case 'pulse'
            ui = a*(on & k < k_start(i) + n_w);
        case 'bipolar'
            ui = a*(on & k < k_start(i) + n_h) - a*(k >= k_start(i) + n_h & k < k_start(i) + n_w);
    end
    u(:, j(i)) = u(:, j(i)) + ui;
end

% the plant's discrete update from x = 0
Ia = zeros(nt, model.n_coils);
Y = zeros(nt, numel(model.y_names));
x = zeros(model.n_states, 1);
stopped_at = NaN;
n_done = nt;
for kk = 1:nt-1
    x = model.A*x + model.B*u(kk,:)';
    Ia(kk+1,:) = (model.x_to_coil*x)';
    Y(kk+1,:) = (model.C_y*x)';
    if max(abs(Ia(kk+1,:))) > seq.i_limit
        stopped_at = kk*seq.dt;
        n_done = kk + 1;
        break
    end
end

[~, name, ext] = fileparts(model_dir);
res.label = [name ext];                           % model names contain a dot
res.model_dir = model_dir;
res.seq = seq;
res.time = k(1:n_done)*seq.dt;
res.supply_names = model.input_names;
res.u = u(1:n_done,:);
res.pulse_start = k_start*seq.dt;
res.coil_names = model.coil_names;
res.Ia = Ia(1:n_done,:);
[res.circuit_names, res.Ic] = circuit_currents(res.Ia, model.coil_names);
res.y_names = model.y_names;
res.Y = Y(1:n_done,:);
res.n_states = model.n_states;
res.growth = max(real(eig(model.A_c)));
res.stopped_at = stopped_at;
end

function seq = check_sequence(seq)
% defaults and checks of the seq struct
defaults = struct('supplies',{{'CS3U','CS2U','CS1U','CS1L','CS2L','CS3L','PF1','PF2', ...
    'PF3','PF4','PF5','VS1','PF6','VS3'}},'waveform','bipolar','amplitude',10, ...
    't_on',0.01,'spacing',2,'width',0.1,'t_end',[],'dt',1e-3,'i_limit',1e6);
if ~isstruct(seq)
    error('supply_sequence_response:input', 'seq must be a struct (help supply_sequence_response).');
end
for f = fieldnames(defaults)'
    if ~isfield(seq, f{1}) || isempty(seq.(f{1}))
        seq.(f{1}) = defaults.(f{1});
    end
end
seq.supplies = reshape(cellstr(seq.supplies), 1, []);
n = numel(seq.supplies);
if isempty(seq.t_end)
    seq.t_end = seq.t_on + n*seq.spacing;
end
seq.waveform = char(seq.waveform);
if ~any(strcmp(seq.waveform, {'bipolar','pulse','step'}))
    error('supply_sequence_response:input', ...
        'waveform must be ''bipolar'', ''pulse'' or ''step'' (got ''%s'').', seq.waveform);
end
amp = double(seq.amplitude(:))';
if isscalar(amp)
    amp = repmat(amp, 1, n);
end
if numel(amp) ~= n || any(~isfinite(amp))
    error('supply_sequence_response:input', ...
        'amplitude must be one value or one per supply (%d supplies, got %d values).', n, numel(amp));
end
seq.amplitude = amp;
for f = {'t_on','spacing','width','t_end','dt','i_limit'}
    v = seq.(f{1});
    if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v)
        error('supply_sequence_response:input', '%s must be a finite number.', f{1});
    end
end
if ~(seq.dt > 0) || seq.t_on < 0 || ~(seq.width > 0) || ~(seq.i_limit > 0)
    error('supply_sequence_response:input', 'Need dt > 0, t_on >= 0, width > 0 and i_limit > 0.');
end
if seq.t_end < seq.t_on + (n-1)*seq.spacing + seq.dt
    error('supply_sequence_response:input', 't_end (%g s) ends before the last supply starts.', seq.t_end);
end
if n > 1 && ~strcmp(seq.waveform,'step') && seq.width > seq.spacing
    error('supply_sequence_response:input', ...
        'width (%g s) must not exceed spacing (%g s): the pulses would overlap.', seq.width, seq.spacing);
end
end
