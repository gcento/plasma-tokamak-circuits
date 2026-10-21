classdef test_plasmaless_linear < matlab.unittest.TestCase
    % Unit tests for the linear plant models (imas_model/configure_linear_model.m). Fixture:
    % linear_models/niceSimpleLin_105073_t84.4176, the circuit form (12 circuit states) of
    % the received NICE linearisation (pulse 105073, confirmed by the user on 2026-10-14,
    % delivered under the label 105078; descriptor with a sign-convention block; raw
    % delivery in linear_models/sources/).
    %
    % Run from the repository root with an IMAS-MATLAB DD-4.x module, reading an
    % IMAS entry first with MATLAB R2024a (see tests/test_plasmaless_imas.m):
    % add_plasmaless_paths; load_plasmaless_md(); runtests('tests/test_plasmaless_linear.m')
    % Nine tests
    % 1) test_build (sizes, power-supply inputs applied to B as they are, circuit
    %    states (CS1U = CS1L, VS3L = -VS3U), sign convention of the descriptor (B
    %    unchanged, Ip row of C unchanged, geometric rows negated), operating point of
    %    the descriptor (IMAS convention), D = 0, stability)
    % 2) test_step_loop (plasmaless_step_ids vs the direct recursion x_{k+1} = A x_k + B u_k;
    %    the 14 power-supply voltages read by index, coil{k}.voltage = supply k)
    % 3) test_operating_point (du = u - op.voltage, output offsets)
    % 4) test_descriptor_errors (broken descriptors raise the documented errors,
    %    incl. the convention block and the ITER circuits; the raw sources are refused;
    %    without convention the file matrices are used as is)
    % 5) test_all_models_circuits (every model of linear_models/, future ones included:
    %    B columns of CS1U and CS1L equal, VS1 = -PF2 -PF3 +PF4 +PF5, signs of the
    %    supply responses, one step with CS1U = CS1L and VS3U = -VS3L)
    % 6) test_observables_round_trip (equilibrium IDS write/read)
    % 7) test_circuit_unchanged (circuit model selections and 4-output step)
    % 8) test_one_pulse_per_generator (input identification: one bipolar pulse per
    %    generator on the continuous file matrices, which currents react)
    % 9) test_absolute_initial_state (linear_models/plasmaless_model: initial state from
    %    the coil currents, error when they violate the circuits)
    properties
        rootDir
        modelDir
        em_coupling
        pf_active
        pf_passive
    end

    methods (TestClassSetup)
        function setup(testCase)
            testCase.rootDir = fileparts(fileparts(mfilename('fullpath')));
            addpath(testCase.rootDir, fullfile(testCase.rootDir,'imas_model'));
            testCase.modelDir = fullfile(testCase.rootDir,'linear_models','niceSimpleLin_105073_t84.4176');
            md_uri = getenv('PLASMALESS_MD_URI');   % empty: data/md_dd4 default
            [testCase.em_coupling, testCase.pf_active, testCase.pf_passive] = ...
                load_plasmaless_md(md_uri);
        end
    end

    methods (Test)

        function test_build(testCase)
            dt = 0.01;
            [model, x0] = configure_linear_model(testCase.modelDir, dt, testCase.pf_active);
            testCase.verifyEqual(model.kind, 'linear');
            testCase.verifyEqual(model.n_states, 16);
            testCase.verifyEqual(model.n_coils, 14);
            testCase.verifyEqual(model.n_loops, 4);
            testCase.verifyEqual(numel(model.y_names), 9);
            testCase.verifySize(model.A, [16 16]);
            testCase.verifySize(model.B, [16 14]);
            testCase.verifySize(model.C_y, [9 16]);
            testCase.verifySize(model.x_to_coil, [14 16]);
            testCase.verifySize(model.x_to_loop, [4 16]);
            testCase.verifySize(x0, [16 1]);
            testCase.verifyEqual(x0, zeros(16,1));
            testCase.verifyEqual(model.coil_names, {'CS3U','CS2U','CS1U','CS1L','CS2L','CS3L', ...
                'PF1','PF2','PF3','PF4','PF5','PF6','VS3U','VS3L'});
            % inputs: the 14 generator voltages (model provider, 2026-10-14)
            testCase.verifyEqual(model.input_names, {'CS3U','CS2U','CS1U','CS1L','CS2L','CS3L', ...
                'PF1','PF2','PF3','PF4','PF5','VS1','PF6','VS3'});
            testCase.verifyEqual(model.state_names, {'CS3U','CS2U','CS1','CS2L','CS3L', ...
                'PF1','PF2','PF3','PF4','PF5','PF6','VS3','PS_1','PS_2','VV_1','VV_2'});
            testCase.verifyEqual(model.loop_names, {'PS_1','PS_2','VV_1','VV_2'});
            % circuits: CS1U and CS1L on the state CS1 (series), VS3U = VS3 and
            % VS3L = -VS3 (antiseries)
            iU = find(strcmp(model.coil_names,'CS1U')); iL = find(strcmp(model.coil_names,'CS1L'));
            testCase.verifyEqual(model.x_to_coil(iU,:), model.x_to_coil(iL,:));
            testCase.verifyEqual(model.x_to_coil(strcmp(model.coil_names,'VS3U'),:), ...
                double(strcmp(model.state_names,'VS3')));
            testCase.verifyEqual(model.x_to_coil(strcmp(model.coil_names,'VS3L'),:), ...
                -double(strcmp(model.state_names,'VS3')));
            testCase.verifyEqual(sum(abs(model.x_to_coil),2), ones(14,1));
            testCase.verifyEqual(sum(model.x_to_loop,2), ones(4,1));
            testCase.verifyEqual(model.D, zeros(9,14));
            testCase.verifyEqual(model.dt, dt);
            testCase.verifyLessThanOrEqual(max(real(eig(model.A_c))), 1e-9);
            % the matrix files, used as they are
            A_npy = read_npy(fullfile(testCase.modelDir,'A.npy'));
            B_npy = read_npy(fullfile(testCase.modelDir,'B.npy'));
            C_npy = read_npy(fullfile(testCase.modelDir,'C.npy'));
            testCase.verifyEqual(A_npy, model.A_c);
            testCase.verifyEqual(model.C_file, C_npy);
            % inputs: the 14 power-supply voltages, applied to B as they are; VS1
            % (column 12) acts as -PF2 -PF3 +PF4 +PF5 (own check, 2026-10-15)
            testCase.verifyEqual(model.B_c, B_npy);
            testCase.verifySize(model.B_c, [16 14]);
            col = @(name) B_npy(:,strcmp(model.input_names,name));
            vs1 = -col('PF2') - col('PF3') + col('PF4') + col('PF5');
            testCase.verifyLessThanOrEqual(norm(col('VS1') - vs1), 1e-12*norm(vs1));
            % sign convention of the descriptor (2026-10-14): state_sign -1,
            % input_sign -1, output_sign Ip -1 and +1 for the geometric outputs.
            % B_conv = (input_sign/state_sign) B = B: B_c is the file B with no sign
            % change; C_conv(k,:) =
            % (state_sign/output_sign(k)) C(k,:): the Ip row is unchanged, the
            % geometric rows (rmag, zmag, rX, zX, gap1..3, Dist2ndSep) are negated.
            testCase.verifyEqual(model.convention.state_sign, -1);
            testCase.verifyEqual(model.convention.input_sign, -1);
            testCase.verifyEqual(model.convention.output_sign, [1;1;1;1;-1;1;1;1;1]);
            testCase.verifyEqual(model.A_file, A_npy);
            testCase.verifyEqual(model.B_file, B_npy);
            testCase.verifyEqual(model.B_c, model.B_file);
            iIp = find(strcmp(model.y_names,'Ip'));
            testCase.verifyEqual(iIp, 5);
            testCase.verifyEqual(model.C_c(iIp,:), C_npy(5,:));
            testCase.verifyEqual(model.C_c(strcmp(model.y_names,'rmag'),:), -C_npy(1,:));
            testCase.verifyEqual(model.C_c, [-C_npy(1:4,:); C_npy(5,:); -C_npy(6:9,:)]);
            testCase.verifyEqual(model.C_y, model.C_c);   % zero-order hold keeps C
            % operating point of the descriptor in the IMAS convention: the negated
            % provider values (model provider, 2026-10-14; frame derived from the DINA
            % 105073 data, 2026-10-14)
            op_cur = -[23248.2 -9514.24 -34331.9 -34331.9 -6293.26 37480.8 21873.6 -29058.6 ...
                -19884.8 -40412.3 -10061.6 19263.3 203.564 -203.564]';
            testCase.verifyEqual(model.op.coil_current, op_cur);
            testCase.verifyEqual(model.op.coil_current_given, true(14,1));
            testCase.verifyEqual(model.op.voltage, zeros(14,1));      % no voltage entry
            testCase.verifyEqual(model.op.voltage_given, false(14,1));
            y_known = ismember(model.y_names, {'rmag','zmag','Ip'})';
            testCase.verifyEqual(model.op.y_known, y_known);
            testCase.verifyEqual(model.op.y(strcmp(model.y_names,'rmag')), 6.29704);
            testCase.verifyEqual(model.op.y(strcmp(model.y_names,'zmag')), 0.596501);
            testCase.verifyEqual(model.op.y(strcmp(model.y_names,'Ip')), -15063200); % provider's 1.50632e7 A, IMAS sign (DINA 105073)
            testCase.verifyEqual(model.source.shot, 105073);
            testCase.verifyEqual(model.op.y(~y_known), zeros(6,1));
        end

        function test_step_loop(testCase)
            dt = 0.01; n_steps = 200;
            model = configure_linear_model(testCase.modelDir, dt, testCase.pf_active);
            % plain deviation model: zero operating point for this test
            model.op.coil_current(:) = 0; model.op.y(:) = 0; model.op.y_known(:) = false;
            TS = timeseries_plasmaless_model(0, n_steps*dt, dt);
            x = zeros(model.n_states,1);
            X = zeros(n_steps+1, model.n_states);
            Ia = zeros(n_steps+1, model.n_coils);
            Ie = zeros(n_steps+1, model.n_loops);
            Yo = zeros(n_steps+1, numel(model.y_names));
            for k = 1:n_steps
                pf_in = test_plasmaless_linear.inputSlice(model, TS.Time(k), TS.Data(k,:));
                [pfa_out, pfp_out, x, y_out] = plasmaless_step_ids(model, pf_in, x);
                X(k+1,:) = x';
                Ia(k+1,:) = cellfun(@(c) c.current.data, pfa_out.coil);
                Ie(k+1,:) = cellfun(@(l) l.current, pfp_out.loop);
                Yo(k+1,:) = y_out.values';
                testCase.assertEqual(numel(pfp_out.loop), 4);
                testCase.assertEqual(y_out.names, model.y_names);
                testCase.assertEqual(y_out.known, false(9,1));
            end
            testCase.verifyEqual(cellfun(@(l) l.name, pfp_out.loop(:)', 'UniformOutput', false), ...
                model.loop_names);
            % direct recursion with the discrete matrices
            Xd = zeros(n_steps+1, model.n_states);
            xd = zeros(model.n_states,1);
            for k = 1:n_steps
                xd = model.A*xd + model.B*TS.Data(k,:)';
                Xd(k+1,:) = xd';
            end
            scale = max(abs(Xd),[],'all');
            testCase.verifyGreaterThan(scale, 0);
            testCase.verifyLessThanOrEqual(max(abs(X - Xd),[],'all'), 1e-12*scale);
            testCase.verifyLessThanOrEqual(max(abs(Ia - Xd*model.x_to_coil'),[],'all'), 1e-12*scale);
            testCase.verifyLessThanOrEqual(max(abs(Ie - Xd*model.x_to_loop'),[],'all'), 1e-12*scale);
            Yd = Xd*model.C_y';
            testCase.verifyLessThanOrEqual(max(abs(Yo - Yd),[],'all'), 1e-12*max(abs(Yd),[],'all'));
            iU = find(strcmp(model.coil_names,'CS1U')); iL = find(strcmp(model.coil_names,'CS1L'));
            testCase.verifyEqual(Ia(:,iU), Ia(:,iL));
            iU = find(strcmp(model.coil_names,'VS3U')); iL = find(strcmp(model.coil_names,'VS3L'));
            testCase.verifyEqual(Ia(:,iU), -Ia(:,iL));
            testCase.verifyTrue(all(isfinite(X(:))));
        end

        function test_operating_point(testCase)
            dt = 0.01;
            model = configure_linear_model(testCase.modelDir, dt, testCase.pf_active);
            model.op.coil_current = (1:14)'*1e3;
            model.op.voltage = (1:14)'*0.5;
            model.op.y = (1:9)'*0.1;
            model.op.y_known = true(9,1);
            x = zeros(model.n_states,1);
            % u = op.voltage: du = 0, the state stays zero and the outputs are op
            pf_in = test_plasmaless_linear.inputSlice(model, 0, model.op.voltage');
            [pfa_out, pfp_out, x1, y_out] = plasmaless_step_ids(model, pf_in, x);
            testCase.verifyEqual(x1, zeros(model.n_states,1));
            testCase.verifyEqual(cellfun(@(c) c.current.data, pfa_out.coil(:)), model.op.coil_current);
            testCase.verifyEqual(cellfun(@(l) l.current, pfp_out.loop(:)), zeros(4,1));
            testCase.verifyEqual(y_out.values, model.op.y);
            testCase.verifyEqual(y_out.known, true(9,1));
            testCase.verifyEqual(cellfun(@(c) c.voltage.data, pfa_out.coil(:)), model.op.voltage);
            % u = op.voltage + v: same state as the plain model driven by v
            v = 10*ones(1,14);
            pf_in = test_plasmaless_linear.inputSlice(model, 0, model.op.voltage' + v);
            [pfa_out, ~, x2, y_out] = plasmaless_step_ids(model, pf_in, x);
            x_ref = model.B*v';
            testCase.verifyEqual(x2, x_ref, 'AbsTol', 1e-12*max(abs(x_ref)));
            testCase.verifyEqual(cellfun(@(c) c.current.data, pfa_out.coil(:)), ...
                model.op.coil_current + model.x_to_coil*x_ref, 'AbsTol', 1e-9);
            testCase.verifyEqual(y_out.values, model.op.y + model.C_y*x_ref, 'AbsTol', 1e-12);
        end

        function test_descriptor_errors(testCase)
            dt = 0.01;
            base = fullfile('/scratch/users', getenv('USER'), 'plasmaless');
            if ~isfolder(base), mkdir(base); end
            d = jsondecode(fileread(fullfile(testCase.modelDir,'linear_model.json')));

            % duplicate input name
            bad = d; bad.inputs{3} = 'CS1L';
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            % input count differs from the columns of B
            bad = d; bad.inputs{end+1} = 'VS1b';
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:size');
            % input_map is no longer accepted: the inputs are the power-supply voltages
            bad = d; bad.input_map = struct('VS3', struct('VS3U', -1, 'VS3L', 1));
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            % inputs: the 14 power supplies in the controller's order (here PF6 and VS1 swapped)
            bad = d; bad.inputs([12 13]) = bad.inputs([13 12]);
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            % missing coil in coil_states
            bad = d; bad.coil_states = rmfield(bad.coil_states,'PF6');
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            % coil mapped to an unknown state
            bad = d; bad.coil_states.PF6 = 'PF7';
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            % unclassified state (VV_2 removed from passive_states)
            bad = d; bad.passive_states = bad.passive_states(1:3);
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            % state both coil and passive
            bad = d; bad.passive_states{end+1} = 'PF6';
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            % coil_states outside the ITER circuits: CS1L on another state, CS1L with
            % the opposite sign of CS1U, VS3L with the same sign as VS3U
            bad = d; bad.coil_states.CS1L = 'CS2L';
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:circuits');
            bad = d; bad.coil_states.CS1L = '-CS1';
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:circuits');
            bad = d; bad.coil_states.VS3L = 'VS3';
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:circuits');
            % the raw sources are not plant models: NICE with VS3 and -VS3 as two
            % states, the plasma-free model with the coil voltages as inputs
            src = fullfile(testCase.rootDir,'linear_models','sources');
            testCase.verifyError(@() configure_linear_model(fullfile(src, ...
                'niceSimpleLin_105073_t84.4176'), dt, testCase.pf_active), 'plasmaless:circuits');
            testCase.verifyError(@() configure_linear_model(fullfile(src,'plasmaless_model'), ...
                dt, testCase.pf_active), 'plasmaless:descriptor');
            % absolute coordinates: circuit states shared by two coils are accepted, an
            % operating point is not
            absolute = rmfield(d,'operating_point'); absolute.coordinates = 'absolute';
            model = testCase.buildFrom(base, absolute, dt);
            testCase.verifyEqual(model.coordinates, 'absolute');
            bad = d; bad.coordinates = 'absolute';
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            % missing required field
            bad = rmfield(d,'outputs');
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            % unknown output name
            bad = d; bad.outputs{1} = 'rgeo';
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:observable');
            % wrong-size matrix
            bad = d; bad.files.A = 'wrongA.txt';
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt, 'wrongA.txt', zeros(15,16)), ...
                'plasmaless:size');
            % missing matrix file
            bad = d; bad.files.B = 'missing.npy';
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:file');
            % operating point with an unknown coil
            bad = d; bad.operating_point.voltage.XX1 = 1;
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            % convention: invalid sign value (not +1/-1), unknown output name in
            % output_sign, invalid output sign, non-object block, unknown field
            bad = d; bad.convention.state_sign = 2;
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            bad = d; bad.convention.input_sign = 'minus';
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            bad = d; bad.convention.output_sign.rgeo = -1;
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            bad = d; bad.convention.output_sign.Ip = 0;
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            bad = d; bad.convention = 'flip';
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            bad = d; bad.convention.current_sign = -1;
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            % no convention block: all signs +1, the file matrices are used unchanged
            none_conv = rmfield(d,'convention');
            model = testCase.buildFrom(base, none_conv, dt);
            B_npy = read_npy(fullfile(testCase.modelDir,'B.npy'));
            C_npy = read_npy(fullfile(testCase.modelDir,'C.npy'));
            testCase.verifyEqual(model.convention.state_sign, 1);
            testCase.verifyEqual(model.convention.input_sign, 1);
            testCase.verifyEqual(model.convention.output_sign, ones(9,1));
            testCase.verifyEqual(model.B_c, B_npy);
            testCase.verifyEqual(model.C_c, C_npy);
            testCase.verifyEqual(model.C_y, C_npy);
            % partial block (only state_sign): B and all C rows negated
            part = none_conv; part.convention = struct('state_sign', -1);
            model = testCase.buildFrom(base, part, dt);
            testCase.verifyEqual(model.convention.input_sign, 1);
            testCase.verifyEqual(model.B_c, -B_npy);
            testCase.verifyEqual(model.C_c, -C_npy);
            % missing descriptor
            testCase.verifyError(@() configure_linear_model(tempname(base), dt, testCase.pf_active), ...
                'plasmaless:descriptor');
            % valid operating point entries are applied and flagged (partial object)
            good = d;
            good.operating_point = struct('voltage', struct('CS3U', 2.5), ...
                'coil_current', struct('PF1', 1e4), 'outputs', struct('Ip', 1.2e7));
            model = testCase.buildFrom(base, good, dt);
            testCase.verifyEqual(model.op.voltage(1), 2.5);
            testCase.verifyEqual(model.op.voltage_given, [true; false(13,1)]);
            testCase.verifyEqual(model.op.coil_current(strcmp(model.coil_names,'PF1')), 1e4);
            testCase.verifyEqual(model.op.coil_current_given, strcmp(model.coil_names,'PF1')');
            testCase.verifyEqual(model.op.y(strcmp(model.y_names,'Ip')), 1.2e7);
            testCase.verifyEqual(model.op.y_known, strcmp(model.y_names,'Ip')');
            % no operating point: zeros, nothing given
            none = rmfield(d,'operating_point');
            model = testCase.buildFrom(base, none, dt);
            testCase.verifyEqual(model.op.coil_current, zeros(14,1));
            testCase.verifyFalse(any(model.op.coil_current_given));
            testCase.verifyFalse(any(model.op.voltage_given));
            testCase.verifyEqual(model.op.y_known, false(9,1));
        end

        function test_all_models_circuits(testCase)
            % every model of linear_models/ (one folder per model; linear_models/sources/
            % holds the raw sources, not plant models) is built as NICE: the 14
            % power-supply voltages in, the ITER circuits as states (checked by
            % configure_linear_model), the coil currents out. New model folders are
            % picked up automatically.
            files = dir(fullfile(testCase.rootDir,'linear_models','*','linear_model.json'));
            testCase.verifyGreaterThanOrEqual(numel(files), 4);
            w = warning('off','plasmaless:unstable');   % vertical models: unstable by design
            restore = onCleanup(@() warning(w));
            for k = 1:numel(files)
                [~, label] = fileparts(files(k).folder);
                model = configure_linear_model(files(k).folder, 0.001, testCase.pf_active);
                % B: CS1U and CS1L drive the same circuit; VS1 = -PF2 -PF3 +PF4 +PF5
                col = @(s) model.B_c(:,strcmp(model.input_names,s));
                testCase.verifyLessThanOrEqual(norm(col('CS1U') - col('CS1L')), ...
                    1e-12*norm(col('CS1L')), [label ': B columns CS1U and CS1L differ']);
                vs1 = -col('PF2') - col('PF3') + col('PF4') + col('PF5');
                testCase.verifyLessThanOrEqual(norm(col('VS1') - vs1), 1e-12*norm(vs1), ...
                    [label ': B column VS1 is not -PF2 -PF3 +PF4 +PF5']);
                % signs (IMAS convention): a positive supply voltage raises the current
                % of its own coil at first; a positive VS3 raises the VS3L current
                dI = model.x_to_coil*model.B_c;     % initial dI/dt per volt, coil x supply
                coil = @(s) strcmp(model.coil_names,s);
                for s = setdiff(model.input_names, {'VS1','VS3'}, 'stable')
                    testCase.verifyGreaterThan(dI(coil(s{1}), strcmp(model.input_names,s{1})), 0, ...
                        sprintf('%s: supply %s does not raise its coil current', label, s{1}));
                end
                j3 = strcmp(model.input_names,'VS3');
                testCase.verifyGreaterThan(dI(coil('VS3L'),j3), 0, [label ': VS3 does not raise the VS3L current']);
                testCase.verifyLessThan(dI(coil('VS3U'),j3), 0, [label ': VS3 does not lower the VS3U current']);
                % one step from the operating point with all supplies driven: series
                % and antiseries coil currents
                u = model.op.voltage' + 100*sin(1:14);
                pfa = plasmaless_step_ids(model, test_plasmaless_linear.inputSlice(model, 0, u), ...
                    zeros(model.n_states,1));
                Ia = cellfun(@(c) c.current.data, pfa.coil(:));
                dIa = Ia - model.op.coil_current;
                tol = 1e-9*max(1, max(abs(Ia)));
                testCase.verifyGreaterThan(max(abs(dIa)), 1e3*tol, [label ': no response']);
                testCase.verifyLessThanOrEqual(abs(dIa(coil('CS1U')) - dIa(coil('CS1L'))), tol, ...
                    [label ': CS1U and CS1L currents differ']);
                testCase.verifyLessThanOrEqual(abs(dIa(coil('VS3U')) + dIa(coil('VS3L'))), tol, ...
                    [label ': VS3U and VS3L currents are not opposite']);
            end
        end

        function test_observables_round_trip(testCase)
            names = {'rmag','zmag','rX','zX','Ip','gap1','gap2','gap3','Dist2ndSep'};
            y.names = names;
            y.values = [6.2; 0.5; 5.1; -3.3; 1.5e7; 0.1; 0.2; 0.3; 0.05];
            y.known = true(9,1);
            eq = ids_init('equilibrium');
            eq.ids_properties.homogeneous_time = 1;
            eq.time = 1;
            eq.time_slice = ids_allocate('equilibrium','time_slice',1);
            eq.time_slice{1}.time = 1;
            eq = plasmaless_observables_to_equilibrium(eq, y);
            s = eq.time_slice{1};
            testCase.verifyEqual(s.global_quantities.magnetic_axis.r, 6.2);
            testCase.verifyEqual(s.global_quantities.ip, 1.5e7);
            testCase.verifyEqual(s.contour_tree.node{1}.critical_type, int32(1));
            testCase.verifyEqual(s.contour_tree.node{1}.r, 5.1);
            testCase.verifyEqual(numel(s.boundary.gap), 3);
            testCase.verifyEqual(s.boundary.gap{3}.name, 'gap3');
            testCase.verifyEqual(s.boundary.gap{3}.value, 0.3);
            [values, known] = plasmaless_observables_from_equilibrium(s, names);
            % Dist2ndSep has no DD4 field: written nowhere, read back unknown
            expected_known = [true(8,1); false];
            testCase.verifyEqual(known, expected_known);
            testCase.verifyEqual(values(expected_known), y.values(expected_known));
            testCase.verifyTrue(isnan(values(end)));
            % partial mask: only the known entries are written
            y2 = y; y2.known = false(9,1); y2.known([2 5 7]) = true;
            eq2 = ids_init('equilibrium');
            eq2.time_slice = ids_allocate('equilibrium','time_slice',1);
            eq2 = plasmaless_observables_to_equilibrium(eq2, y2);
            [v2, k2] = plasmaless_observables_from_equilibrium(eq2.time_slice{1}, names);
            testCase.verifyEqual(k2, y2.known);
            testCase.verifyEqual(v2(k2), y.values(k2));
            testCase.verifyEqual(numel(eq2.time_slice{1}.boundary.gap), 2);
            % the fresh slice has nothing known
            [~, k0] = plasmaless_observables_from_equilibrium( ...
                ids_allocate('equilibrium','time_slice',1), names);
            testCase.verifyEqual(k0, false(9,1));
            testCase.verifyError(@() plasmaless_observables_from_equilibrium(s, {'rgeo'}), ...
                'plasmaless:observable');
        end

        function test_circuit_unchanged(testCase)
            dt = 0.01;
            model = configure_plasmaless_model(testCase.em_coupling, testCase.pf_active, ...
                testCase.pf_passive, dt);
            n = model.n_coils + model.n_loops;
            testCase.verifyEqual(model.kind, 'circuit');
            testCase.verifyEqual(model.n_states, n);
            testCase.verifyEqual(model.state_names, [model.coil_names model.loop_names]);
            testCase.verifyEqual(model.x_to_coil, [eye(model.n_coils) zeros(model.n_coils,model.n_loops)]);
            testCase.verifyEqual(model.x_to_loop, [zeros(model.n_loops,model.n_coils) eye(model.n_loops)]);
            testCase.verifyEmpty(model.y_names);
            testCase.verifySize(model.C_y, [0 n]);
            testCase.verifyEqual(model.op.voltage, zeros(model.n_coils,1));
            x = zeros(n,1);
            TS = timeseries_plasmaless_model(0, 20*dt, dt);
            for k = 1:20
                pf_in = test_plasmaless_linear.inputSlice(model, TS.Time(k), TS.Data(k,:));
                [pfa3, pfp3, x3] = plasmaless_step_ids(model, pf_in, x);
                [pfa4, pfp4, x4, y_out] = plasmaless_step_ids(model, pf_in, x);
                testCase.assertEqual(x4, x3);
                testCase.assertEqual(cellfun(@(c) c.current.data, pfa4.coil), ...
                    cellfun(@(c) c.current.data, pfa3.coil));
                testCase.assertEqual(cellfun(@(l) l.current, pfp4.loop), ...
                    cellfun(@(l) l.current, pfp3.loop));
                % exactly the state, as before the generalisation
                testCase.assertEqual(cellfun(@(c) c.current.data, pfa4.coil(:)), x4(1:model.n_coils));
                testCase.assertEqual(cellfun(@(l) l.current, pfp4.loop(:)), x4(model.n_coils+1:end));
                testCase.assertEqual(x4, model.A*x + model.B*TS.Data(k,:)');
                testCase.assertEmpty(y_out.values);
                x = x4;
            end
            testCase.verifyEqual(pfa4.ids_properties.comment, 'plasmaless coil-vessel model: coil currents');
        end

        function test_one_pulse_per_generator(testCase)
            % input identification in the provider's generator space: the
            % continuous file matrices (A_file, B_file, provider frame; the
            % checks are on relative patterns, so no sign convention is needed),
            % c2d at dt = 1e-3, x0 = 0; per generator +10 V for 0.5 s, -10 V
            % for 0.5 s, then 0.5 s of zero
            dt = 1e-3;
            model = configure_linear_model(testCase.modelDir, 0.01, testCase.pf_active);
            [M, F, X] = test_plasmaless_linear.pulsePerGenerator(model.A_file, model.B_file, dt);
            gen = model.input_names; st = model.state_names;
            % 1) each CS/PF generator drives its own state hardest among the 11
            %    CS/PF coil states, with a positive initial response
            own = {'CS3U','CS3U'; 'CS2U','CS2U'; 'CS1U','CS1'; 'CS1L','CS1'; 'CS2L','CS2L'; ...
                'CS3L','CS3L'; 'PF1','PF1'; 'PF2','PF2'; 'PF3','PF3'; 'PF4','PF4'; ...
                'PF5','PF5'; 'PF6','PF6'};
            for r = 1:size(own,1)
                j = find(strcmp(gen, own{r,1})); i = find(strcmp(st, own{r,2}));
                [~, imax] = max(M(1:11,j));
                testCase.verifyEqual(st{imax}, own{r,2}, sprintf( ...
                    'generator %s: largest CS/PF response on state %s, expected %s', ...
                    own{r,1}, st{imax}, own{r,2}));
                testCase.verifyGreaterThan(F(i,j), 0, sprintf( ...
                    'generator %s: initial response of state %s is not positive (%g A)', ...
                    own{r,1}, own{r,2}, F(i,j)));
                testCase.verifyGreaterThan(model.B_file(i,j), 0, sprintf( ...
                    'generator %s: B_file(%s,%s) is not positive', own{r,1}, own{r,2}, own{r,1}));
            end
            % 2) CS1U and CS1L feed the same state: identical trajectories
            jU = find(strcmp(gen,'CS1U')); jL = find(strcmp(gen,'CS1L'));
            scale = max(abs(X{jU}),[],'all');
            testCase.verifyLessThanOrEqual(max(abs(X{jU} - X{jL}),[],'all'), 1e-12*scale, ...
                'generators CS1U and CS1L: state trajectories differ');
            % 3) VS1 (star point in series with PF2-PF5): measured max |deviation|
            %    (2026-10-14) PF2 13.13, PF3 2.905, PF4 4.504, PF5 3.608 A, initial
            %    signs - - + +; CS coils at most 0.869 A (CS2U), PF1 0.592, PF6 0.473 A
            %    (VS3 about 83 A, not part of this check)
            jV = find(strcmp(gen,'VS1'));
            iStar = find(ismember(st, {'PF2','PF3','PF4','PF5'}));
            iOther = find(ismember(st, {'CS3U','CS2U','CS1','CS2L','CS3L','PF1','PF6'}));
            for i = iStar
                for o = iOther
                    testCase.verifyGreaterThan(M(i,jV), M(o,jV), sprintf( ...
                        'generator VS1: state %s (%g A) does not exceed state %s (%g A)', ...
                        st{i}, M(i,jV), st{o}, M(o,jV)));
                end
            end
            % 4) VS3: the circuit state VS3 (the VS3U current, VS3L = -VS3) responds
            %    most among the 12 circuit states, negative at first (B(12,14) < 0)
            j3 = find(strcmp(gen,'VS3'));
            [~, imax] = max(M(1:12,j3));
            testCase.verifyEqual(st{imax}, 'VS3', sprintf( ...
                'generator VS3: largest circuit-state response on %s, expected VS3', st{imax}));
            testCase.verifyLessThan(F(12,j3), 0, 'generator VS3: initial response of state VS3 is not negative');
        end

        function test_absolute_initial_state(testCase)
            % linear_models/plasmaless_model (absolute coordinates, circuit states): the
            % initial state built from coil currents that follow the circuits gives them
            % back, with passive states 0; currents that violate the circuits are refused
            model = configure_linear_model(fullfile(testCase.rootDir,'linear_models', ...
                'plasmaless_model'), 0.01, testCase.pf_active);
            testCase.verifyEqual(model.coordinates, 'absolute');
            coil = @(s) strcmp(model.coil_names,s);
            i0 = 1e3*(1:14)';
            i0(coil('CS1L')) = i0(coil('CS1U'));        % series
            i0(coil('VS3L')) = -i0(coil('VS3U'));       % antiseries
            x0 = linear_initial_state(model, i0);
            testCase.verifySize(x0, [model.n_states 1]);
            testCase.verifyEqual(model.x_to_coil*x0, i0, 'AbsTol', 1e-12*max(abs(i0)));
            testCase.verifyEqual(model.x_to_loop*x0, zeros(model.n_loops,1));
            bad = i0; bad(coil('CS1L')) = bad(coil('CS1L')) + 100;
            testCase.verifyError(@() linear_initial_state(model, bad), 'plasmaless:circuits');
            bad = i0; bad(coil('VS3L')) = -bad(coil('VS3L'));
            testCase.verifyError(@() linear_initial_state(model, bad), 'plasmaless:circuits');
        end

    end

    methods (Access = private)
        function model = buildFrom(testCase, base, desc, dt, extra_name, extra_matrix)
            % build from a temporary copy of the model directory with descriptor desc
            dir_tmp = tempname(base);
            copyfile(testCase.modelDir, dir_tmp);
            fid = fopen(fullfile(dir_tmp,'linear_model.json'),'w');
            % Inf/NaN written as Infinity/NaN (jsondecode reads them back)
            fprintf(fid,'%s', jsonencode(desc,'ConvertInfAndNaN',false));
            fclose(fid);
            if nargin > 4
                writematrix(extra_matrix, fullfile(dir_tmp,extra_name), 'Delimiter',' ', 'FileType','text');
            end
            cleanup = onCleanup(@() rmdir(dir_tmp,'s'));
            model = configure_linear_model(dir_tmp, dt, testCase.pf_active);
        end
    end

    methods (Static, Access = private)
        function [M, F, X] = pulsePerGenerator(A, B, dt)
            % one bipolar pulse per generator column of B on the continuous
            % (A, B), discretised with c2d (zero-order hold) at dt, from x0 = 0:
            % +10 V for 0.5 s, -10 V for 0.5 s, 0.5 s of zero. M(i,j) = max |x_i|
            % for generator j, F(:,j) = x after the first step, X{j} = trajectory
            n = size(A,1); m = size(B,2);
            sysd = c2d(ss(A, B, eye(n), zeros(n,m)), dt);
            n1 = round(0.5/dt);
            u = [10*ones(1,n1), -10*ones(1,n1), zeros(1,n1)];
            M = zeros(n,m); F = zeros(n,m); X = cell(1,m);
            for j = 1:m
                x = zeros(n,1); Xj = zeros(numel(u), n);
                for k = 1:numel(u)
                    x = sysd.A*x + sysd.B(:,j)*u(k);
                    Xj(k,:) = x';
                end
                M(:,j) = max(abs(Xj),[],1)';
                F(:,j) = Xj(1,:)';
                X{j} = Xj;
            end
        end

        function pf_in = inputSlice(model, t, u)
            % one-slice pf_active with voltages u at time t in coil{k}.voltage (for a
            % linear model u(k) is power supply k; coils named as the model coils)
            pf_in = ids_init('pf_active');
            pf_in.ids_properties.homogeneous_time = 1;
            pf_in.time = t;
            pf_in.coil = ids_allocate('pf_active','coil',model.n_coils);
            for ii = 1:model.n_coils
                pf_in.coil{ii}.name = model.coil_names{ii};
                pf_in.coil{ii}.voltage.data = u(ii);
            end
        end
    end
end
