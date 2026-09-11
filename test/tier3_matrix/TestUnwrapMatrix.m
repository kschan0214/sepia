%% TestUnwrapMatrix - Tier 3 real-dataset phase-unwrapping method regression matrix
%
% Mirrors TestQSMMatrix.m: runs the one-stop SEPIA pipeline on a real
% dataset once per phase-unwrapping method in methodUnwrapName
% (sepia_universal_variables.m), holding background field removal
% ('VSHARP') and QSM dipole inversion ('TKD', both toolbox-free) fixed so
% the unwrap method is the only isolated variable. '3D best path' is
% deliberately excluded (see sepiatest.unwrap_method_toolbox_key.m).
%
classdef TestUnwrapMatrix < matlab.unittest.TestCase

    properties (TestParameter)
        unwrapMethod = get_unwrap_methods();
    end

    properties
        Dataset
        Toolboxes
    end

    methods (TestClassSetup)
        function setupPathAndData(testCase)
            thisFile   = mfilename('fullpath');
            testRoot   = fileparts(fileparts(thisFile));
            SEPIA_HOME = fileparts(testRoot);

            addpath(SEPIA_HOME);
            sepia_addpath;
            addpath(testRoot);

            testCase.Dataset   = sepiatest.get_real_dataset();
            testCase.Toolboxes = sepiatest.discover_toolboxes();
            addpath(testRoot);

            testCase.assumeTrue(~isempty(testCase.Dataset.inputDir) && ~isempty(testCase.Dataset.maskFile), ...
                ['Tier 3 real-dataset tests skipped: set SEPIA_TEST_REAL_DATA_DIR/SEPIA_TEST_REAL_DATA_MASK ', ...
                 '(or test/config/real_dataset.json) to point at a real dataset. See test/README.md.']);
        end
    end

    methods (Test)
        function testUnwrapMethodMatchesReference(testCase, unwrapMethod)
            import matlab.unittest.fixtures.TemporaryFolderFixture

            addpath(fileparts(fileparts(mfilename('fullpath'))));
            toolboxKey = sepiatest.unwrap_method_toolbox_key(unwrapMethod);
            sepiatest.assume_toolbox_available(testCase, testCase.Toolboxes, toolboxKey, unwrapMethod);

            fixture = testCase.applyFixture(TemporaryFolderFixture);
            outputPrefix = fullfile(fixture.Folder, 'sepia');

            algorParam = struct();
            algorParam.general.isInvert = false;
            algorParam.general.isBET    = false;
            algorParam.unwrap.unwrapMethod   = unwrapMethod;
            algorParam.unwrap.echoCombMethod = 'Optimum weights';
            algorParam.bfr.method = 'VSHARP';
            algorParam.qsm.method = 'TKD';

            sepiaIO(testCase.Dataset.inputDir, outputPrefix, testCase.Dataset.maskFile, algorParam);
            addpath(fileparts(fileparts(mfilename('fullpath'))));

            suffix = '.nii.gz';
            slug   = matlab.lang.makeValidName(unwrapMethod);

            sepiatest.assert_output_contract(testCase, outputPrefix, suffix, ...
                { {'fieldmap', true}, {'Chimap', true} }, testCase.Dataset.maskFile);

            maskImg   = load_nii_img_only(testCase.Dataset.maskFile) > 0;
            fieldImg  = load_nii_img_only([outputPrefix '_fieldmap' suffix]);
            actual    = sepiatest.mask_stats(fieldImg, maskImg);

            refFile = fullfile(sepiatest.test_root(), 'references', 'tier3', sprintf('unwrap_%s.mat', slug));
            testCase.assumeTrue(isfile(refFile), ...
                sprintf(['No Tier-3 reference saved yet for unwrap method "%s". ', ...
                         'Run regenerate_reference(''confirm'',true,''tier'',''tier3'') to create it after a manual review.'], unwrapMethod));

            ref = load(refFile);
            tol = sepiatest.tolerance_for_method('unwrap', unwrapMethod);
            sepiatest.compare_to_reference(testCase, actual, ref.stats, tol, sprintf('Tier2:Unwrap:%s', unwrapMethod));
        end
    end

end

function methods = get_unwrap_methods()
thisFile   = mfilename('fullpath');
testRoot   = fileparts(fileparts(thisFile));
SEPIA_HOME = fileparts(testRoot);
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
methods = methodUnwrapName(:)';
end
