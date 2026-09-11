%% TestQSMMatrixPhantom - Tier 2 synthetic-phantom QSM method regression matrix
%
% Runs the one-stop SEPIA pipeline on a small, deterministically-generated
% synthetic phantom (test/phantom/generate_synthetic_phantom.m) once per
% QSM dipole-inversion method in methodQSMName (sepia_universal_variables.m)
% AND per matrix size (the same 5 odd/even sizes TestOddMatrixSize.m uses
% in Tier 1), holding phase unwrapping ('None') and background field
% removal ('VSHARP', the built-in variant) fixed so the QSM method is the
% only isolated variable.
%
% Unlike Tier 1 (test/tier1_smoke/), this exercises EVERY method, including
% toolbox-dependent ones - the phantom is small (32x32x24 by default, vs. a
% real dataset's ~176x256x144), so even iterative solvers finish quickly.
% This also closes a gap TestOddMatrixSize.m documents: some pad/crop paths
% (e.g. Laplacian-based unwrap) only trigger for toolbox-dependent methods,
% unreachable in Tier 1's toolbox-free suite.
%
% Each row is skipped (not failed) when the method's required external
% toolbox isn't installed, or it has no distributed model checkpoints
% anywhere in this environment (see qsm_method_toolbox_key.m). Only the
% baseline matrix size ([32 32 24]) is compared numerically against a saved
% reference; the odd-size cases only assert output shape/no-NaN (see
% assert_output_contract.m) - they're about pad/crop correctness, not
% algorithm-accuracy regression, so no reference is needed per size.
%
classdef TestQSMMatrixPhantom < matlab.unittest.TestCase

    properties (TestParameter)
        qsmMethod  = get_qsm_methods();
        matrixSize = {[32 32 24], [31 32 24], [32 31 24], [32 32 25], [31 31 25]};
    end

    properties
        Toolboxes
    end

    methods (TestClassSetup)
        function setupPathAndToolboxes(testCase)
            thisFile   = mfilename('fullpath');          % .../test/tier2_matrix/TestQSMMatrixPhantom
            testRoot   = fileparts(fileparts(thisFile));  % .../test
            SEPIA_HOME = fileparts(testRoot);             % repo root

            addpath(SEPIA_HOME);
            sepia_addpath; % auto-creates SpecifyToolboxesDirectory.m if missing
            addpath(testRoot);
            addpath(fullfile(testRoot, 'phantom'));

            testCase.Toolboxes = sepiatest.discover_toolboxes();
            addpath(testRoot); addpath(fullfile(testRoot, 'phantom')); % discover_toolboxes strips test/ off the path
        end
    end

    methods (Test)
        function testQSMMethodMatchesReference(testCase, qsmMethod, matrixSize)
            import matlab.unittest.fixtures.TemporaryFolderFixture

            addpath(fileparts(fileparts(mfilename('fullpath')))); % defensive: ensure test/ is still on path
            toolboxKey = sepiatest.qsm_method_toolbox_key(qsmMethod);
            sepiatest.assume_toolbox_available(testCase, testCase.Toolboxes, toolboxKey, qsmMethod);

            fixture = testCase.applyFixture(TemporaryFolderFixture);
            phantomPaths = generate_synthetic_phantom(fullfile(fixture.Folder, 'phantom'), matrixSize);

            outFixture   = testCase.applyFixture(TemporaryFolderFixture);
            outputPrefix = fullfile(outFixture.Folder, 'sepia');

            algorParam = struct();
            algorParam.general.isInvert = false;
            algorParam.general.isBET    = false;
            algorParam.unwrap.unwrapMethod   = 'None';
            algorParam.unwrap.echoCombMethod = 'Optimum weights';
            algorParam.bfr.method = 'VSHARP';
            algorParam.qsm.method = qsmMethod;
            if strcmp(qsmMethod, 'NDI')
                algorParam.qsm.isGPU = false; % force CPU-only for determinism
            end

            sepiaIO(phantomPaths.input, outputPrefix, phantomPaths.mask, algorParam);

            % sepiaIO strips test/, this file's own tier2_matrix/ folder,
            % and test/phantom/ back off the path via sepia_addpath - see
            % test/tier1_smoke/TestOddMatrixSize.m for why all three matter.
            addpath(fileparts(fileparts(mfilename('fullpath'))));
            addpath(fileparts(mfilename('fullpath')));
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'phantom'));

            suffix = '.nii.gz';
            slug   = matlab.lang.makeValidName(qsmMethod);

            sepiatest.assert_output_contract(testCase, outputPrefix, suffix, ...
                { {'Chimap', true}, {'localfield', true}, {'fieldmap', true} }, ...
                phantomPaths.mask);

            if ~isequal(matrixSize, [32 32 24])
                return % odd-size cases: shape/no-NaN check above is sufficient, no numeric reference per size
            end

            chiFile = [outputPrefix '_Chimap' suffix];
            maskImg = load_nii_img_only(phantomPaths.mask) > 0;
            chiImg  = load_nii_img_only(chiFile);
            actual  = sepiatest.mask_stats(chiImg, maskImg);

            refFile = fullfile(sepiatest.test_root(), 'references', 'tier2', sprintf('qsm_%s.mat', slug));
            testCase.assumeTrue(isfile(refFile), ...
                sprintf(['No Tier-2 reference saved yet for QSM method "%s". ', ...
                         'Run regenerate_reference(''confirm'',true,''tier'',''tier2'') to create it after a manual review.'], qsmMethod));

            ref = load(refFile);
            tol = sepiatest.tolerance_for_method('qsm', qsmMethod);
            sepiatest.compare_to_reference(testCase, actual, ref.stats, tol, sprintf('Tier2:QSM:%s', qsmMethod));
        end
    end

end

%% pull the method list from SEPIA itself rather than hardcoding it here
function methods = get_qsm_methods()

thisFile   = mfilename('fullpath');           % .../test/tier2_matrix/TestQSMMatrixPhantom
testRoot   = fileparts(fileparts(thisFile));  % .../test
SEPIA_HOME = fileparts(testRoot);             % repo root

if exist('sepia_universal_variables','file') ~= 2
    % bare addpath(SEPIA_HOME) is not enough - sepia_universal_variables
    % itself needs configuration/ (and other subfolders) on path too,
    % which only the real sepia_addpath sets up. This matters because
    % TestSuite.fromFolder evaluates TestParameter defaults (i.e. calls
    % this function) at suite-CONSTRUCTION time, before any
    % TestClassSetup method has run.
    addpath(SEPIA_HOME);
    sepia_addpath;
end
sepia_universal_variables;
methods = methodQSMName(:)';

end
