classdef test_generators < matlab.unittest.TestCase
    % The generators of generators/ reproduce the stored files:
    % 1) test_plasmaless_export: export_plasmaless_model (data/md_dd4) =
    %    generators/plasmaless/source;
    % 2) test_plasmaless_transform: transform_plasmaless_model of the stored source =
    %    linear_models/plasmaless_model;
    % 3) test_nice_transform: transform_nice_lin_model of generators/nice/sources/<name> =
    %    linear_models/<name>, for niceLin and niceSimpleLin.
    % Same descriptor fields and matrices within 1e-12 (utils/same_model.m).
    % The vertical model needs the PDS equilibrium: utils/make_vertical_circuit_model.m
    % with compare_with checks it. The folders are written under
    % /scratch/users/$USER/linear_model and removed at the end.
    % Run from the repository root with an IMAS-MATLAB DD-4.x module (README.md):
    % add_linear_model_paths; load_machine_description(); runtests('tests/test_generators.m')
    properties
        rootDir
        outBase
    end

    methods (TestClassSetup)
        function setup(testCase)
            testCase.rootDir = fileparts(fileparts(mfilename('fullpath')));
            gen = fullfile(testCase.rootDir,'generators');
            addpath(fullfile(testCase.rootDir,'plant'), fullfile(testCase.rootDir,'utils'), ...
                fullfile(gen,'common'), fullfile(gen,'plasmaless'), fullfile(gen,'nice'));
            base = fullfile('/scratch/users', getenv('USER'), 'linear_model');
            if ~isfolder(base), mkdir(base); end
            testCase.outBase = tempname(base);
            mkdir(testCase.outBase);
            testCase.addTeardown(@rmdir, testCase.outBase, 's');
        end
    end

    methods (Test)
        function test_plasmaless_export(testCase)
            out = fullfile(testCase.outBase,'plasmaless_source');
            export_plasmaless_model(out);
            testCase.verifySameModel(out, fullfile(testCase.rootDir,'generators','plasmaless','source'));
        end

        function test_plasmaless_transform(testCase)
            out = fullfile(testCase.outBase,'plasmaless_model');
            transform_plasmaless_model(fullfile(testCase.rootDir,'generators','plasmaless','source'), out);
            testCase.verifySameModel(out, fullfile(testCase.rootDir,'linear_models','plasmaless_model'));
        end

        function test_nice_transform(testCase)
            for name = {'niceLin_105073_t84.4176', 'niceSimpleLin_105073_t84.4176'}
                out = fullfile(testCase.outBase, name{1});
                transform_nice_lin_model(fullfile(testCase.rootDir,'generators','nice','sources',name{1}), out);
                testCase.verifySameModel(out, fullfile(testCase.rootDir,'linear_models',name{1}));
            end
        end
    end

    methods (Access = private)
        function verifySameModel(testCase, dir_new, dir_ref)
            [same, msg] = same_model(dir_new, dir_ref);
            testCase.verifyTrue(same, msg);
        end
    end
end
