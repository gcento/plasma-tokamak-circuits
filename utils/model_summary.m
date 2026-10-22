function model = model_summary(model_dir, pf_active_md)
% Load a model folder with the plant's loader (configure_linear_model, dt = 1 ms) and
% print what the plant sees: states (circuits + passive), coordinates, outputs, growth
% rate and the circuit checks (CS1U = CS1L, VS3U = -VS3L, B column of VS1 = -PF2 -PF3
% +PF4 +PF5). The loader itself refuses a folder that breaks the contract of
% linear_models/README.md.
%   model = model_summary(model_dir)
%   model = model_summary(model_dir, pf_active_md)   % machine-description pf_active

if nargin < 2 || isempty(pf_active_md)
    [~, pf_active_md] = load_machine_description();
end
w = warning('off','linear_model:unstable');        % vertical models: unstable by design
restore = onCleanup(@() warning(w));
model = configure_linear_model(char(model_dir), 1e-3, pf_active_md);
row = @(name) model.x_to_coil(strcmp(model.coil_names,name),:);
col = @(name) model.B_c(:,strcmp(model.input_names,name));
vs1 = -col('PF2') - col('PF3') + col('PF4') + col('PF5');
fprintf(['%s\n  %d states (%d circuits + %d passive), coordinates %s, outputs {%s}\n', ...
    '  growth rate %.4g 1/s; CS1U = CS1L: %d, VS3U = -VS3L: %d, ', ...
    'VS1 = -PF2 -PF3 +PF4 +PF5: %d\n'], char(model_dir), model.n_states, ...
    model.n_states - model.n_loops, model.n_loops, model.coordinates, ...
    strjoin(model.y_names(:)',','), max(real(eig(model.A_c))), ...
    isequal(row('CS1U'),row('CS1L')), isequal(row('VS3U'),-row('VS3L')), ...
    norm(col('VS1') - vs1) <= 1e-12*norm(vs1));
end
