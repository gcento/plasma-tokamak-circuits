% run_supply_sequence.m
% Open-loop test of the plant models: voltages on the power supplies, one supply after
% the other (utils/supply_sequence_response.m), and the response of the 12 circuit
% currents in time, one model per colour. With one supply it is the single-channel test.
% Optionally overlays the same test run through the MUSCLE3 actor
% (muscle3/run_actor_standalone.sh): the two must coincide.
% Output: one line per model, one figure (supply voltage + 12 circuit currents), one
% figure with zmag for the models with a plasma, and optionally a .mat file.
% Every curve is a deviation from the operating point of its model. In an open-loop
% unstable model (niceLin, vertical_circuit) the vertical mode grows after a few tenths
% of a second: keep the run short (t_end) or use the stable models.
% Edit the parameters, then run it in MATLAB on a compute node (IMAS-MATLAB module loaded).

%% parameters
close all
clear
clc

models    = {'plasmaless_model','vertical_circuit_105073_t84.3176','niceLin_105073_t84.4176'};
                            % folder names of linear_models/ or absolute paths, or 'all'
supplies  = {'VS3'};        % in pulse order, any of CS3U CS2U CS1U CS1L CS2L CS3L PF1 PF2
                            % PF3 PF4 PF5 VS1 PF6 VS3
waveform  = 'bipolar';      % 'bipolar' (+A then -A) | 'pulse' (+A) | 'step' (+A until the end)
amplitude = -100;           % [V] one value, or one per supply
t_on      = 0.01;           % [s] start of the first supply (>= dt to compare with MUSCLE3)
spacing   = 0.3;            % [s] between two supplies: long enough for the response to die
                            %     out (plasmaless_model: 0.5 % of the peak left after 1 s)
width     = 0.2;            % [s] duration of pulse and bipolar
t_end     = [];             % [s] end of the run ([] = t_on + number of supplies * spacing)
dt        = 1e-3;           % [s] time step (plant.dt)
i_limit   = 1e6;            % [A] stop a model when a coil current exceeds it
md_uri    = 'imas:hdf5?path=/home/ITER/centong/Desktop/PDSs/IMAS-PDS/local_installs/plasmaless-tokamak-circuits/data/md_dd4';             % machine description ('' = data/md_dd4 of this repository)
do_plot   = true;
results_file    = './VS3_niceLin_vs_md.mat';   % .mat to save the results in ('' = not saved)
muscle3_results = '';       % channel_response.mat of run_actor_standalone.sh ('' = none):
                            % same model, one supply, same waveform

%% run
here = fileparts(mfilename('fullpath'));
repo_root = fileparts(here);
addpath(fullfile(repo_root,'plant'), here);
if ischar(models) && strcmp(models,'all')
    files = dir(fullfile(repo_root,'linear_models','*','linear_model.json'));
    model_dirs = {files.folder};
else
    model_dirs = cellstr(models);
    for m = 1:numel(model_dirs)
        if ~startsWith(model_dirs{m}, filesep)
            model_dirs{m} = fullfile(repo_root,'linear_models',model_dirs{m});
        end
    end
end
[~, pf_active_md] = load_machine_description(md_uri);
seq = struct('supplies',{cellstr(supplies)},'waveform',waveform,'amplitude',amplitude, ...
    't_on',t_on,'spacing',spacing,'width',width,'t_end',t_end,'dt',dt,'i_limit',i_limit);
results = cell(1, numel(model_dirs));
for m = 1:numel(model_dirs)
    results{m} = supply_sequence_response(model_dirs{m}, seq, pf_active_md);
    r = results{m};
    fprintf('%s: %d states, growth %.3g 1/s', r.label, r.n_states, r.growth);
    if ~isnan(r.stopped_at)
        fprintf(', stopped at t = %.3g s (coil current > i_limit)', r.stopped_at);
    end
    fprintf('\n');
end
if ~isempty(muscle3_results)
    results{end+1} = read_muscle3_results(muscle3_results, results{1}.seq);
    compare_with_matlab(results{end}, results(1:end-1));
end
if do_plot
    plot_supply_sequence(results);
end
if ~isempty(results_file)
    save(results_file, 'results', '-v7');
    fprintf('Results saved to %s\n', results_file);
end

%% local functions
function r = read_muscle3_results(file, seq)
% the results file of muscle3/run_actor_standalone.sh in the form of
% supply_sequence_response (single supply)
s = load(file);
same = numel(seq.supplies) == 1 && strcmp(char(s.channel), seq.supplies{1}) && ...
    strcmp(char(s.waveform), seq.waveform) && s.amplitude == seq.amplitude(1) && ...
    abs(s.t_on - seq.t_on) < 1e-9 && abs(s.width - seq.width) < 1e-9 && abs(s.dt - seq.dt) < 1e-12;
if ~same
    warning('run_supply_sequence:muscle3', ['%s was run with other inputs (%s %s %g V, ', ...
        't_on %g s, width %g s, dt %g s).'], file, char(s.channel), char(s.waveform), ...
        s.amplitude, s.t_on, s.width, s.dt);
end
[~, name, ext] = fileparts(char(s.model_dir));
r.label = ['MUSCLE3: ' name ext];
r.seq = seq;
r.time = s.time(:);
supply_names = {'CS3U','CS2U','CS1U','CS1L','CS2L','CS3L','PF1','PF2','PF3','PF4','PF5', ...
    'VS1','PF6','VS3'};
r.u = zeros(numel(r.time), numel(supply_names));
r.u(:, strcmp(supply_names, char(s.channel))) = s.u(:);
r.pulse_start = seq.t_on;
r.coil_names = cellstr(s.coil_names);
r.Ia = s.Ia;
[r.circuit_names, r.Ic] = circuit_currents(r.Ia, r.coil_names);
r.y_names = cellstr(s.y_names);
r.Y = s.Y;
r.stopped_at = NaN;
r.line_style = '--';
end

function compare_with_matlab(r3, results)
% one line: largest difference of the circuit currents from the MATLAB result of the
% same model, relative to its largest circuit current
model = extractAfter(r3.label, 'MUSCLE3: ');
k = find(strcmp(cellfun(@(r) r.label, results, 'UniformOutput', false), model), 1);
if isempty(k)
    fprintf('MUSCLE3 run of %s: model not among models, not compared\n', model);
    return
end
n = min(numel(r3.time), numel(results{k}.time));
ref = results{k}.Ic(1:n,:);
d = max(abs(r3.Ic(1:n,:) - ref), [], 'all')/max(max(abs(ref), [], 'all'), realmin);
fprintf('MUSCLE3 run of %s: largest difference from MATLAB %.2g (relative)\n', model, d);
end
