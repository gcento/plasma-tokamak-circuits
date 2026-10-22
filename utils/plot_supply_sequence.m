function plot_supply_sequence(results)
% Plot the results of utils/supply_sequence_response.m, one colour per model (a result
% with the field line_style, e.g. '--' for MUSCLE3, uses it):
%   1) the power-supply voltage on top (each supply named at its start) and the 12
%      circuit currents in time below;
%   2) zmag of the models with a plasma.
% The time axes are linked: zooming one tile zooms all.
%   plot_supply_sequence(results)               % results: cell array of res structs

if isstruct(results)
    results = num2cell(results);
end
n = numel(results);
colors = lines(max(n,1));
labels = cellfun(@(r) strrep(r.label, '_', '\_'), results, 'UniformOutput', false);
for m = 1:n
    if ~isnan(results{m}.stopped_at)
        labels{m} = sprintf('%s (stopped at %.3g s)', labels{m}, results{m}.stopped_at);
    end
end
[~, longest] = max(cellfun(@(r) numel(r.time), results));
rv = results{longest};                            % the inputs are the same for every model
v = sum(rv.u, 2);                                 % one supply at a time
v_max = max(max(abs(v)), eps);

% 1) supply voltage and circuit currents
figure('Name', 'supply sequence');
tl = tiledlayout(5, 3, 'TileSpacing', 'compact');
title(tl, sprintf('%s on the power supplies, one every %g s', rv.seq.waveform, rv.seq.spacing));
ax = nexttile(tl, [1 3]);
stairs(rv.time, v, 'LineWidth', 1.5);
for i = 1:numel(rv.seq.supplies)
    text(rv.pulse_start(i), 1.2*v_max, rv.seq.supplies{i}, 'FontSize', 8);
end
ylim([-1.5 1.5]*v_max);
grid on
title('Power-supply voltage');
ylabel('V [V]');
for c = 1:numel(rv.circuit_names)
    ax(end+1) = nexttile(tl); %#ok<AGROW>
    hold on
    for m = 1:n
        plot(results{m}.time, results{m}.Ic(:,c), line_style(results{m}), 'Color', colors(m,:), 'LineWidth', 1);
    end
    grid on
    title(rv.circuit_names{c});
    ylabel('\DeltaI [A]');
end
xlabel(tl, 't [s]');
lg = legend(ax(end), labels);
lg.Layout.Tile = 'east';

% 2) zmag of the models with a plasma
has_z = cellfun(@(r) any(strcmp(r.y_names, 'zmag')) && ...
    any(isfinite(r.Y(2:end, strcmp(r.y_names, 'zmag')))), results);
if any(has_z)
    figure('Name', 'zmag');
    ax(end+1) = axes;
    hold on
    for m = find(has_z)
        r = results{m};
        plot(r.time, 1e3*r.Y(:, strcmp(r.y_names, 'zmag')), line_style(r), 'Color', colors(m,:), 'LineWidth', 1);
    end
    grid on
    title('Plasma vertical position (magnetic axis)');
    xlabel('t [s]');
    ylabel('\Deltaz [mm]');
    legend(labels(has_z));
end
linkaxes(ax, 'x');
xlim(ax(1), [0 rv.time(end)]);
end

function s = line_style(r)
% '-' unless the result sets line_style
s = '-';
if isfield(r, 'line_style')
    s = r.line_style;
end
end
