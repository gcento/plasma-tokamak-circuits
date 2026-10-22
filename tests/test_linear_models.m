classdef test_linear_models < matlab.unittest.TestCase
    % Tests of the plant (plant/) and of the models of linear_models/. Fixture:
    % linear_models/niceSimpleLin_105073_t84.4176 (12 circuit states + 4 passive, with a
    % sign-convention block).
    % Run from the repository root with an IMAS-MATLAB DD-4.x module (README.md):
    %   add_linear_model_paths; load_machine_description(); runtests('tests/test_linear_models.m')
    % (load_machine_description first: with MATLAB R2024a an IMAS entry must be read before the
    % first runtests call, HDF5 1.8 of the test framework vs 1.14 of the access layer).
    %  1) test_build: sizes, names, circuits, sign convention and operating point
    %  2) test_step_loop: plant_step_ids = the recursion x(k+1) = A x(k) + B u(k)
    %  3) test_operating_point: du = u - op.voltage, outputs = op + deviation
    %  4) test_descriptor_errors: broken descriptors raise the documented errors
    %  5) test_descriptor_variants: convention block, operating point, coordinates
    %  6) test_all_models_circuits: every folder of linear_models/ follows the ITER circuits
    %  7) test_observables_round_trip: model outputs written to and read from equilibrium
    %  8) test_absolute_initial_state: plasmaless_model starts from the coil currents
    %  9) test_supply_sequence: utils/supply_sequence_response.m
    properties
        rootDir
        modelDir
        tmpBase
        pf_active
    end

    methods (TestClassSetup)
        function setup(testCase)
            testCase.rootDir = fileparts(fileparts(mfilename('fullpath')));
            addpath(fullfile(testCase.rootDir,'plant'), fullfile(testCase.rootDir,'utils'));
            testCase.modelDir = fullfile(testCase.rootDir,'linear_models','niceSimpleLin_105073_t84.4176');
            testCase.tmpBase = fullfile('/scratch/users', getenv('USER'), 'linear_model');
            if ~isfolder(testCase.tmpBase), mkdir(testCase.tmpBase); end
            md_uri = getenv('LINEAR_MODEL_MD_URI');   % empty: data/md_dd4
            [~, testCase.pf_active] = load_machine_description(md_uri);
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
            testCase.verifyEqual(model.input_names, {'CS3U','CS2U','CS1U','CS1L','CS2L','CS3L', ...
                'PF1','PF2','PF3','PF4','PF5','VS1','PF6','VS3'});
            testCase.verifyEqual(model.state_names, {'CS3U','CS2U','CS1','CS2L','CS3L', ...
                'PF1','PF2','PF3','PF4','PF5','PF6','VS3','PS_1','PS_2','VV_1','VV_2'});
            testCase.verifyEqual(model.loop_names, {'PS_1','PS_2','VV_1','VV_2'});
            % circuits: CS1U and CS1L on the state CS1, VS3U = VS3 and VS3L = -VS3
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
            % A and B as in the files; VS1 acts as -PF2 -PF3 +PF4 +PF5
            A_npy = read_npy(fullfile(testCase.modelDir,'A.npy'));
            B_npy = read_npy(fullfile(testCase.modelDir,'B.npy'));
            C_npy = read_npy(fullfile(testCase.modelDir,'C.npy'));
            testCase.verifyEqual(A_npy, model.A_c);
            testCase.verifyEqual(model.C_file, C_npy);
            testCase.verifyEqual(model.B_c, B_npy);
            testCase.verifySize(model.B_c, [16 14]);
            col = @(name) B_npy(:,strcmp(model.input_names,name));
            vs1 = -col('PF2') - col('PF3') + col('PF4') + col('PF5');
            testCase.verifyLessThanOrEqual(norm(col('VS1') - vs1), 1e-12*norm(vs1));
            % convention of the descriptor: state_sign = input_sign = -1 (B unchanged),
            % output_sign -1 for Ip (Ip row of C unchanged, geometric rows negated)
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
            % operating point of the descriptor in the IMAS convention (the provider's
            % values negated; Ip -15.0632 MA as in DINA 105073)
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
            testCase.verifyEqual(model.op.y(strcmp(model.y_names,'Ip')), -15063200);
            testCase.verifyEqual(model.source.shot, 105073);
            testCase.verifyEqual(model.op.y(~y_known), zeros(6,1));
        end

        function test_step_loop(testCase)
            dt = 0.01; n_steps = 200;
            model = configure_linear_model(testCase.modelDir, dt, testCase.pf_active);
            % plain deviation model: zero operating point for this test
            model.op.coil_current(:) = 0; model.op.y(:) = 0; model.op.y_known(:) = false;
            [tt, U] = test_linear_models.pulseTrain(0, n_steps*dt, dt);
            x = zeros(model.n_states,1);
            X = zeros(n_steps+1, model.n_states);
            Ia = zeros(n_steps+1, model.n_coils);
            Ie = zeros(n_steps+1, model.n_loops);
            Yo = zeros(n_steps+1, numel(model.y_names));
            for k = 1:n_steps
                pf_in = test_linear_models.inputSlice(model, tt(k), U(k,:));
                [pfa_out, pfp_out, x, y_out] = plant_step_ids(model, pf_in, x);
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
                xd = model.A*xd + model.B*U(k,:)';
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
            pf_in = test_linear_models.inputSlice(model, 0, model.op.voltage');
            [pfa_out, pfp_out, x1, y_out] = plant_step_ids(model, pf_in, x);
            testCase.verifyEqual(x1, zeros(model.n_states,1));
            testCase.verifyEqual(cellfun(@(c) c.current.data, pfa_out.coil(:)), model.op.coil_current);
            testCase.verifyEqual(cellfun(@(l) l.current, pfp_out.loop(:)), zeros(4,1));
            testCase.verifyEqual(y_out.values, model.op.y);
            testCase.verifyEqual(y_out.known, true(9,1));
            testCase.verifyEqual(cellfun(@(c) c.voltage.data, pfa_out.coil(:)), model.op.voltage);
            % u = op.voltage + v: same state as the plain model driven by v
            v = 10*ones(1,14);
            pf_in = test_linear_models.inputSlice(model, 0, model.op.voltage' + v);
            [pfa_out, ~, x2, y_out] = plant_step_ids(model, pf_in, x);
            x_ref = model.B*v';
            testCase.verifyEqual(x2, x_ref, 'AbsTol', 1e-12*max(abs(x_ref)));
            testCase.verifyEqual(cellfun(@(c) c.current.data, pfa_out.coil(:)), ...
                model.op.coil_current + model.x_to_coil*x_ref, 'AbsTol', 1e-9);
            testCase.verifyEqual(y_out.values, model.op.y + model.C_y*x_ref, 'AbsTol', 1e-12);
        end

        function test_descriptor_errors(testCase)
            d = jsondecode(fileread(fullfile(testCase.modelDir,'linear_model.json')));
            % change to the fixture descriptor, expected error
            cases = {
                @(d) setfield(d,'inputs',strrep(d.inputs,'CS1U','CS1L')),      'linear_model:descriptor' % duplicate input
                @(d) setfield(d,'inputs',[d.inputs; {'VS1b'}]),                'linear_model:size'       % 15 inputs, 14 columns of B
                @(d) setfield(d,'inputs',d.inputs([1:11 13 12 14])),           'linear_model:descriptor' % not the controller's order
                @(d) setfield(d,'input_map',struct('VS3U',-1,'VS3L',1)),       'linear_model:descriptor' % input_map no longer used
                @(d) setfield(d,'coil_states',rmfield(d.coil_states,'PF6')),   'linear_model:descriptor' % coil missing
                @(d) setfield(d,'coil_states','PF6','PF7'),                    'linear_model:descriptor' % unknown state
                @(d) setfield(d,'passive_states',d.passive_states(1:3)),       'linear_model:descriptor' % VV_2 neither coil nor passive
                @(d) setfield(d,'passive_states',[d.passive_states; {'PF6'}]), 'linear_model:descriptor' % PF6 both coil and passive
                @(d) setfield(d,'coil_states','CS1L','CS2L'),                  'linear_model:circuits'   % CS1L not with CS1U
                @(d) setfield(d,'coil_states','CS1L','-CS1'),                  'linear_model:circuits'   % CS1L opposite to CS1U
                @(d) setfield(d,'coil_states','VS3L','VS3'),                   'linear_model:circuits'   % VS3L not in antiseries
                @(d) setfield(d,'coordinates','absolute'),                     'linear_model:descriptor' % absolute with an operating point
                @(d) rmfield(d,'outputs'),                                     'linear_model:descriptor' % required field missing
                @(d) setfield(d,'outputs',strrep(d.outputs,'rmag','rgeo')),    'linear_model:observable' % output not in the vocabulary
                @(d) setfield(d,'files','B','missing.npy'),                    'linear_model:file'       % matrix file missing
                @(d) setfield(d,'operating_point','voltage','XX1',1),          'linear_model:descriptor' % unknown power supply
                @(d) setfield(d,'convention','state_sign',2),                  'linear_model:descriptor' % sign not +1/-1
                @(d) setfield(d,'convention','input_sign','minus'),            'linear_model:descriptor' % sign not a number
                @(d) setfield(d,'convention','output_sign','rgeo',-1),         'linear_model:descriptor' % not an output
                @(d) setfield(d,'convention','output_sign','Ip',0),            'linear_model:descriptor' % sign not +1/-1
                @(d) setfield(d,'convention','flip'),                          'linear_model:descriptor' % not an object
                @(d) setfield(d,'convention','current_sign',-1),               'linear_model:descriptor' % unknown field
                };
            for c = 1:size(cases,1)
                bad = cases{c,1}(d);
                testCase.verifyError(@() testCase.buildFrom(bad), cases{c,2}, func2str(cases{c,1}));
            end
            % matrix of the wrong size, missing descriptor
            bad = d; bad.files.A = 'wrongA.txt';
            testCase.verifyError(@() testCase.buildFrom(bad, 'wrongA.txt', zeros(15,16)), ...
                'linear_model:size');
            testCase.verifyError(@() configure_linear_model(tempname(testCase.tmpBase), 0.01, ...
                testCase.pf_active), 'linear_model:descriptor');
            % the raw sources of generators/ are not plant models: NICE with VS3 and -VS3
            % as two states, the plasma-free model with the coil voltages as inputs
            gen = fullfile(testCase.rootDir,'generators');
            testCase.verifyError(@() configure_linear_model(fullfile(gen,'nice','sources', ...
                'niceSimpleLin_105073_t84.4176'), 0.01, testCase.pf_active), 'linear_model:circuits');
            testCase.verifyError(@() configure_linear_model(fullfile(gen,'plasmaless','source'), ...
                0.01, testCase.pf_active), 'linear_model:descriptor');
        end

        function test_descriptor_variants(testCase)
            % valid changes to the fixture descriptor
            d = jsondecode(fileread(fullfile(testCase.modelDir,'linear_model.json')));
            B_npy = read_npy(fullfile(testCase.modelDir,'B.npy'));
            C_npy = read_npy(fullfile(testCase.modelDir,'C.npy'));
            % no convention block: every sign +1, the file matrices used as they are
            model = testCase.buildFrom(rmfield(d,'convention'));
            testCase.verifyEqual(model.convention.state_sign, 1);
            testCase.verifyEqual(model.convention.input_sign, 1);
            testCase.verifyEqual(model.convention.output_sign, ones(9,1));
            testCase.verifyEqual(model.B_c, B_npy);
            testCase.verifyEqual(model.C_c, C_npy);
            testCase.verifyEqual(model.C_y, C_npy);
            % only state_sign = -1: B and every row of C negated
            part = d; part.convention = struct('state_sign', -1);
            model = testCase.buildFrom(part);
            testCase.verifyEqual(model.convention.input_sign, 1);
            testCase.verifyEqual(model.B_c, -B_npy);
            testCase.verifyEqual(model.C_c, -C_npy);
            % partial operating point: the given entries are applied and flagged
            good = d;
            good.operating_point = struct('voltage', struct('CS3U', 2.5), ...
                'coil_current', struct('PF1', 1e4), 'outputs', struct('Ip', 1.2e7));
            model = testCase.buildFrom(good);
            testCase.verifyEqual(model.op.voltage(1), 2.5);
            testCase.verifyEqual(model.op.voltage_given, [true; false(13,1)]);
            testCase.verifyEqual(model.op.coil_current(strcmp(model.coil_names,'PF1')), 1e4);
            testCase.verifyEqual(model.op.coil_current_given, strcmp(model.coil_names,'PF1')');
            testCase.verifyEqual(model.op.y(strcmp(model.y_names,'Ip')), 1.2e7);
            testCase.verifyEqual(model.op.y_known, strcmp(model.y_names,'Ip')');
            % no operating point: zeros, nothing given
            none = rmfield(d,'operating_point');
            model = testCase.buildFrom(none);
            testCase.verifyEqual(model.op.coil_current, zeros(14,1));
            testCase.verifyFalse(any(model.op.coil_current_given));
            testCase.verifyFalse(any(model.op.voltage_given));
            testCase.verifyEqual(model.op.y_known, false(9,1));
            % absolute coordinates, without operating point
            absolute = none; absolute.coordinates = 'absolute';
            model = testCase.buildFrom(absolute);
            testCase.verifyEqual(model.coordinates, 'absolute');
        end

        function test_all_models_circuits(testCase)
            % every folder of linear_models/, new ones included, is built as NICE: the 14
            % power-supply voltages in, the ITER circuits as states, the coil currents out
            files = dir(fullfile(testCase.rootDir,'linear_models','*','linear_model.json'));
            testCase.verifyGreaterThanOrEqual(numel(files), 4);
            w = warning('off','linear_model:unstable');   % vertical models: unstable by design
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
                pfa = plant_step_ids(model, test_linear_models.inputSlice(model, 0, u), ...
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
            eq = plant_observables_to_equilibrium(eq, y);
            s = eq.time_slice{1};
            testCase.verifyEqual(s.global_quantities.magnetic_axis.r, 6.2);
            testCase.verifyEqual(s.global_quantities.ip, 1.5e7);
            testCase.verifyEqual(s.contour_tree.node{1}.critical_type, int32(1));
            testCase.verifyEqual(s.contour_tree.node{1}.r, 5.1);
            testCase.verifyEqual(numel(s.boundary.gap), 3);
            testCase.verifyEqual(s.boundary.gap{3}.name, 'gap3');
            testCase.verifyEqual(s.boundary.gap{3}.value, 0.3);
            [values, known] = plant_observables_from_equilibrium(s, names);
            % Dist2ndSep has no DD4 field: written nowhere, read back unknown
            expected_known = [true(8,1); false];
            testCase.verifyEqual(known, expected_known);
            testCase.verifyEqual(values(expected_known), y.values(expected_known));
            testCase.verifyTrue(isnan(values(end)));
            % partial mask: only the known entries are written
            y2 = y; y2.known = false(9,1); y2.known([2 5 7]) = true;
            eq2 = ids_init('equilibrium');
            eq2.time_slice = ids_allocate('equilibrium','time_slice',1);
            eq2 = plant_observables_to_equilibrium(eq2, y2);
            [v2, k2] = plant_observables_from_equilibrium(eq2.time_slice{1}, names);
            testCase.verifyEqual(k2, y2.known);
            testCase.verifyEqual(v2(k2), y.values(k2));
            testCase.verifyEqual(numel(eq2.time_slice{1}.boundary.gap), 2);
            % the fresh slice has nothing known
            [~, k0] = plant_observables_from_equilibrium( ...
                ids_allocate('equilibrium','time_slice',1), names);
            testCase.verifyEqual(k0, false(9,1));
            testCase.verifyError(@() plant_observables_from_equilibrium(s, {'rgeo'}), ...
                'linear_model:observable');
        end

        function test_absolute_initial_state(testCase)
            % linear_models/plasmaless_model: the initial state built from coil currents
            % that follow the circuits gives them back, with passive states 0; currents
            % that violate the circuits are refused
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
            testCase.verifyError(@() linear_initial_state(model, bad), 'linear_model:circuits');
            bad = i0; bad(coil('VS3L')) = -bad(coil('VS3L'));
            testCase.verifyError(@() linear_initial_state(model, bad), 'linear_model:circuits');
        end

        function test_supply_sequence(testCase)
            % utils/supply_sequence_response.m on the fixture: three supplies one after
            % the other, bipolar pulses of 50 samples every 100 samples
            seq = struct('supplies',{{'PF2','VS1','CS1U'}},'waveform','bipolar', ...
                'amplitude',[10 -20 5],'t_on',0.01,'spacing',0.1,'width',0.05,'dt',1e-3);
            res = supply_sequence_response(testCase.modelDir, seq, testCase.pf_active);
            nt = round((seq.t_on + 3*seq.spacing)/seq.dt) + 1;     % default t_end
            k = (0:nt-1)';
            testCase.verifyEqual(res.time, k*seq.dt);
            testCase.verifySize(res.Ia, [nt 14]);
            testCase.verifySize(res.Ic, [nt 12]);
            testCase.verifyEqual(res.pulse_start, [0.01; 0.11; 0.21], 'AbsTol', 1e-12);
            testCase.verifyTrue(isnan(res.stopped_at));
            % voltages on the sample index, as muscle3/voltage_driver.py: +A for 25
            % samples from k0, then -A for 25 samples; every other supply at 0 V
            bipolar = @(a, k0) a*(k >= k0 & k < k0 + 25) - a*(k >= k0 + 25 & k < k0 + 50);
            u = @(s) res.u(:, strcmp(res.supply_names, s));
            testCase.verifyEqual([u('PF2') u('VS1') u('CS1U')], ...
                [bipolar(10,10) bipolar(-20,110) bipolar(5,210)]);
            testCase.verifyEqual(nnz(res.u), 150);
            % the trajectory of the actor's step function (plant_step_ids) with a
            % zero operating point
            model = configure_linear_model(testCase.modelDir, seq.dt, testCase.pf_active);
            model.op.coil_current(:) = 0; model.op.y(:) = 0;
            x = zeros(model.n_states,1);
            Ia = zeros(nt,14); Y = zeros(nt,numel(model.y_names));
            for kk = 1:nt-1
                pf_in = test_linear_models.inputSlice(model, res.time(kk), res.u(kk,:));
                [pfa, ~, x, y_out] = plant_step_ids(model, pf_in, x);
                Ia(kk+1,:) = cellfun(@(c) c.current.data, pfa.coil);
                Y(kk+1,:) = y_out.values';
            end
            scale = max(abs(Ia),[],'all');
            testCase.verifyGreaterThan(scale, 0);
            testCase.verifyLessThanOrEqual(max(abs(res.Ia - Ia),[],'all'), 1e-12*scale);
            testCase.verifyLessThanOrEqual(max(abs(res.Y - Y),[],'all'), 1e-12*max(abs(Y),[],'all'));
            % circuit currents as the controller measures them: here CS1 = CS1U, VS3 = VS3U
            coil = @(c) res.Ia(:, strcmp(res.coil_names, c));
            testCase.verifyEqual(res.Ic(:, strcmp(res.circuit_names,'CS1')), coil('CS1U'));
            testCase.verifyEqual(res.Ic(:, strcmp(res.circuit_names,'VS3')), coil('VS3U'));
            % superposition: the run is the sum of one-supply runs started at the same times
            Ia = zeros(nt,14);
            for i = 1:3
                one = seq;
                one.supplies = seq.supplies(i);
                one.amplitude = seq.amplitude(i);
                one.t_on = res.pulse_start(i);
                one.t_end = res.time(end);
                r = supply_sequence_response(testCase.modelDir, one, testCase.pf_active);
                Ia = Ia + r.Ia;
            end
            testCase.verifyLessThanOrEqual(max(abs(res.Ia - Ia),[],'all'), 1e-12*scale);
            % step and pulse on one supply
            one = struct('supplies',{{'PF2'}},'amplitude',10,'t_on',0.01,'width',0.05, ...
                't_end',0.2,'dt',1e-3);
            k1 = (0:200)';
            one.waveform = 'step';
            r = supply_sequence_response(testCase.modelDir, one, testCase.pf_active);
            testCase.verifyEqual(r.u(:, strcmp(r.supply_names,'PF2')), 10*(k1 >= 10));
            one.waveform = 'pulse';
            r = supply_sequence_response(testCase.modelDir, one, testCase.pf_active);
            testCase.verifyEqual(r.u(:, strcmp(r.supply_names,'PF2')), 10*(k1 >= 10 & k1 < 60));
            % i_limit stops the run
            stop = seq; stop.i_limit = 1e-3*scale;
            r = supply_sequence_response(testCase.modelDir, stop, testCase.pf_active);
            testCase.verifyFalse(isnan(r.stopped_at));
            testCase.verifyLessThan(numel(r.time), nt);
            % errors: unknown supply, unknown waveform, overlapping pulses, t_end before
            % the last supply starts
            bad = {setfield(seq,'supplies',{'PF2','PF7','CS1U'}), setfield(seq,'waveform','sine'), ...
                setfield(seq,'width',0.2), setfield(seq,'t_end',0.15)};
            for b = 1:numel(bad)
                testCase.verifyError(@() supply_sequence_response(testCase.modelDir, bad{b}, ...
                    testCase.pf_active), 'supply_sequence_response:input');
            end
        end

    end

    methods (Access = private)
        function model = buildFrom(testCase, desc, matrix_file, matrix)
            % build from a temporary copy of the fixture with descriptor desc (and an
            % extra matrix file)
            dir_tmp = tempname(testCase.tmpBase);
            copyfile(testCase.modelDir, dir_tmp);
            cleanup = onCleanup(@() rmdir(dir_tmp,'s'));
            fid = fopen(fullfile(dir_tmp,'linear_model.json'),'w');
            % Inf/NaN written as Infinity/NaN (jsondecode reads them back)
            fprintf(fid,'%s', jsonencode(desc,'ConvertInfAndNaN',false));
            fclose(fid);
            if nargin > 2
                writematrix(matrix, fullfile(dir_tmp,matrix_file), 'Delimiter',' ', 'FileType','text');
            end
            model = configure_linear_model(dir_tmp, 0.01, testCase.pf_active);
        end
    end

    methods (Static, Access = private)
        function [t, U] = pulseTrain(tstart, tend, dt)
            % one bipolar pulse per power supply in turn (14 columns): +10 V for 0.5 s,
            % -10 V for 0.5 s, 1 s of zero between the supplies, on t = (tstart:dt:tend)'
            t = (tstart:dt:tend)';
            U = zeros(numel(t), 14);
            for k = 1:14
                t_start = (k-1)*2;
                U(t >= t_start & t < t_start + 0.5, k) = 10;
                U(t >= t_start + 0.5 & t < t_start + 1, k) = -10;
            end
        end

        function pf_in = inputSlice(model, t, u)
            % one-slice pf_active with the power-supply voltages u at time t,
            % coil{k}.voltage = supply k (coils named as the model coils)
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
