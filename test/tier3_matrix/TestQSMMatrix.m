%% TestQSMMatrix - Tier 3 real-dataset QSM method regression matrix
%
% Runs the one-stop SEPIA pipeline on a real dataset (path supplied via
% sepiatest.get_real_dataset(), never hardcoded) once per QSM
% dipole-inversion method in methodQSMName (sepia_universal_variables.m),
% holding phase unwrapping ('None') and background field removal
% ('VSHARP', the built-in variant) fixed so the QSM method is the only
% isolated variable. Each row is skipped (not failed) when:
%   - no real dataset is configured on this machine, or
%   - the method's required external toolbox isn't installed, or
%   - the method has no distributed model checkpoints anywhere in this
%     environment (deep-learning add-ons - see qsm_method_toolbox_key.m)
%
classdef TestQSMMatrix < matlab.unittest.TestCase

    properties (TestParameter)
        qsmMethod = get_qsm_methods();
    end

    properties
        Dataset
        Toolboxes
    end

    methods (TestClassSetup)
        function setupPathAndData(testCase)
            thisFile   = mfilename('fullpath');          % .../test/tier3_matrix/TestQSMMatrix
            testRoot   = fileparts(fileparts(thisFile));  % .../test
            SEPIA_HOME = fileparts(testRoot);             % repo root

            addpath(SEPIA_HOME);
            sepia_addpath; % auto-creates SpecifyToolboxesDirectory.m if missing
            addpath(testRoot);

            testCase.Dataset   = sepiatest.get_real_dataset();
            testCase.Toolboxes = sepiatest.discover_toolboxes();
            addpath(testRoot); % discover_toolboxes calls sepia_addpath internally, which strips test/ off the path

            testCase.assumeTrue(~isempty(testCase.Dataset.inputDir) && ~isempty(testCase.Dataset.maskFile), ...
                ['Tier 3 real-dataset tests skipped: set SEPIA_TEST_REAL_DATA_DIR/SEPIA_TEST_REAL_DATA_MASK ', ...
                 '(or test/config/real_dataset.json) to point at a real dataset. See test/README.md.']);
        end
    end

    methods (Test)
        function testQSMMethodMatchesReference(testCase, qsmMethod)
            import matlab.unittest.fixtures.TemporaryFolderFixture

            addpath(fileparts(fileparts(mfilename('fullpath')))); % defensive: ensure test/ is still on path
            toolboxKey = sepiatest.qsm_method_toolbox_key(qsmMethod);
            sepiatest.assume_toolbox_available(testCase, testCase.Toolboxes, toolboxKey, qsmMethod);

            fixture = testCase.applyFixture(TemporaryFolderFixture);
            outputPrefix = fullfile(fixture.Folder, 'sepia');

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

            sepiaIO(testCase.Dataset.inputDir, outputPrefix, testCase.Dataset.maskFile, algorParam);

            % sepiaIO/SepiaIOWrapper calls sepia_addpath internally, which
            % strips test/ back off the path - re-add before using
            % sepiatest.* below.
            addpath(fileparts(fileparts(mfilename('fullpath'))));

            suffix = '.nii.gz';
            slug   = matlab.lang.makeValidName(qsmMethod);

            sepiatest.assert_output_contract(testCase, outputPrefix, suffix, ...
                { {'Chimap', true}, {'localfield', true}, {'fieldmap', true} }, ...
                testCase.Dataset.maskFile);

            chiFile = [outputPrefix '_Chimap' suffix];
            maskImg = load_nii_img_only(testCase.Dataset.maskFile) > 0;
            chiImg  = load_nii_img_only(chiFile);
            actual  = sepiatest.mask_stats(chiImg, maskImg);

            refFile = fullfile(sepiatest.test_root(), 'references', 'tier3', sprintf('qsm_%s.mat', slug));
            testCase.assumeTrue(isfile(refFile), ...
                sprintf(['No Tier-3 reference saved yet for QSM method "%s". ', ...
                         'Run regenerate_reference(''confirm'',true,''tier'',''tier3'') to create it after a manual review.'], qsmMethod));

            ref = load(refFile);
            tol = sepiatest.tolerance_for_method('qsm', qsmMethod);
            sepiatest.compare_to_reference(testCase, actual, ref.stats, tol, sprintf('Tier2:QSM:%s', qsmMethod));
        end
    end

end

%% pull the method list from SEPIA itself rather than hardcoding it here
function methods = get_qsm_methods()

thisFile   = mfilename('fullpath');           % .../test/tier3_matrix/TestQSMMatrix
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
