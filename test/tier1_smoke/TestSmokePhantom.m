%% TestSmokePhantom - Tier 1 toolbox-free smoke test
%
% Runs the full one-stop SEPIA pipeline (sepiaIO -> SepiaIOWrapper) on a
% small deterministic synthetic phantom, using only methods that require
% no external toolbox (unwrap='None', bfr='VSHARP' built-in, qsm one of
% {TKD, Closed-form solution, iLSQR}). Confirms the pipeline runs without
% error, produces the expected output files with the expected JSON
% sidecars, and matches saved reference summary statistics within a tight
% tolerance (all three methods are deterministic/closed-form).
%
classdef TestSmokePhantom < matlab.unittest.TestCase

    properties (TestParameter)
        qsmMethod = {'TKD', 'Closed-form solution', 'iLSQR'};
    end

    properties
        PhantomPaths
    end

    methods (TestClassSetup)
        function addSepiaToPath(testCase)
            % Compute paths via mfilename rather than sepiatest.* here,
            % since sepia_addpath (below) does `rmpath(genpath(SEPIA_HOME))`
            % internally, which would strip test/ back off the path if it
            % had already been added by the caller (e.g. run_tier1.m) -
            % so test/ must be (re-)added *after* sepia_addpath runs, not
            % relied upon before it.
            thisFile   = mfilename('fullpath');          % .../test/tier1_smoke/TestSmokePhantom
            testRoot   = fileparts(fileparts(thisFile)); % .../test
            SEPIA_HOME = fileparts(testRoot);            % repo root

            addpath(SEPIA_HOME);
            sepia_addpath; % sets up SpecifyToolboxesDirectory.m (auto-created if missing) and core subfolders

            addpath(testRoot);
            addpath(fullfile(testRoot, 'phantom'));
        end

        function setupPhantom(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            fixture = testCase.applyFixture(TemporaryFolderFixture);
            testCase.PhantomPaths = generate_synthetic_phantom(fullfile(fixture.Folder, 'phantom'));
        end
    end

    methods (Test)
        function testOneStopBuiltinPipeline(testCase, qsmMethod)
            import matlab.unittest.fixtures.TemporaryFolderFixture

            fixture = testCase.applyFixture(TemporaryFolderFixture);
            outDir  = fixture.Folder;
            outputPrefix = fullfile(outDir, 'sepia');

            algorParam = struct();
            algorParam.general.isInvert = false;
            algorParam.general.isBET    = false;
            algorParam.unwrap.unwrapMethod   = 'None';
            algorParam.unwrap.echoCombMethod = 'Optimum weights';
            algorParam.bfr.method  = 'VSHARP';
            algorParam.qsm.method  = qsmMethod;

            sepiaIO(testCase.PhantomPaths.input, outputPrefix, testCase.PhantomPaths.mask, algorParam);

            % sepiaIO/SepiaIOWrapper calls sepia_addpath internally, which
            % does rmpath(genpath(SEPIA_HOME)) and so strips test/ AND
            % this file's own tier1_smoke/ folder back off the path -
            % both must be re-added: test/ before using any sepiatest.*
            % helper below, and tier1_smoke/ or this class itself becomes
            % undispatchable for the next parameterized test case.
            addpath(fileparts(fileparts(mfilename('fullpath'))));
            addpath(fileparts(mfilename('fullpath')));

            suffix = '.nii.gz';
            sepiatest.assert_output_contract(testCase, outputPrefix, suffix, ...
                { {'Chimap', true}, {'localfield', true}, {'fieldmap', true} }, ...
                testCase.PhantomPaths.mask);

            chiFile  = [outputPrefix '_Chimap' suffix];
            maskImg  = load_nii_img_only(testCase.PhantomPaths.mask) > 0;
            chiImg   = load_nii_img_only(chiFile);
            actual   = sepiatest.mask_stats(chiImg, maskImg);

            refFile = fullfile(sepiatest.test_root(), 'references', 'tier1', ...
                sprintf('smoke_%s.mat', matlab.lang.makeValidName(qsmMethod)));
            testCase.assumeTrue(isfile(refFile), ...
                sprintf(['No Tier-1 reference saved yet for QSM method "%s". ', ...
                         'Run test/tools/regenerate_reference.m to create it after a manual review.'], qsmMethod));

            ref = load(refFile);
            tol = sepiatest.tolerance_for_method('qsm', qsmMethod);
            sepiatest.compare_to_reference(testCase, actual, ref.stats, tol, sprintf('Tier1:QSM:%s', qsmMethod));
        end
    end

end
