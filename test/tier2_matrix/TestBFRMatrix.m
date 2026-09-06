%% TestBFRMatrix - Tier 2 real-dataset background-field-removal method regression matrix
%
% Mirrors TestQSMMatrix.m: runs the one-stop SEPIA pipeline on a real
% dataset once per background-field-removal method in methodBFRName
% (sepia_universal_variables.m), holding phase unwrapping ('None') and
% QSM dipole inversion ('TKD', toolbox-free) fixed so the BFR method is
% the only isolated variable.
%
classdef TestBFRMatrix < matlab.unittest.TestCase

    properties (TestParameter)
        bfrMethod = get_bfr_methods();
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
                ['Tier 2 real-dataset tests skipped: set SEPIA_TEST_REAL_DATA_DIR/SEPIA_TEST_REAL_DATA_MASK ', ...
                 '(or test/config/real_dataset.json) to point at a real dataset. See test/README.md.']);
        end
    end

    methods (Test)
        function testBFRMethodMatchesReference(testCase, bfrMethod)
            import matlab.unittest.fixtures.TemporaryFolderFixture

            addpath(fileparts(fileparts(mfilename('fullpath'))));
            toolboxKey = sepiatest.bfr_method_toolbox_key(bfrMethod);
            sepiatest.assume_toolbox_available(testCase, testCase.Toolboxes, toolboxKey, bfrMethod);

            fixture = testCase.applyFixture(TemporaryFolderFixture);
            outputPrefix = fullfile(fixture.Folder, 'sepia');

            algorParam = struct();
            algorParam.general.isInvert = false;
            algorParam.general.isBET    = false;
            algorParam.unwrap.unwrapMethod   = 'None';
            algorParam.unwrap.echoCombMethod = 'Optimum weights';
            algorParam.bfr.method = bfrMethod;
            algorParam.qsm.method = 'TKD';

            sepiaIO(testCase.Dataset.inputDir, outputPrefix, testCase.Dataset.maskFile, algorParam);
            addpath(fileparts(fileparts(mfilename('fullpath'))));

            suffix = '.nii.gz';
            slug   = matlab.lang.makeValidName(bfrMethod);

            sepiatest.assert_output_contract(testCase, outputPrefix, suffix, ...
                { {'localfield', true}, {'Chimap', true} }, testCase.Dataset.maskFile);

            maskImg = load_nii_img_only(testCase.Dataset.maskFile) > 0;
            localFieldImg = load_nii_img_only([outputPrefix '_localfield' suffix]);
            actual = sepiatest.mask_stats(localFieldImg, maskImg);

            refFile = fullfile(sepiatest.test_root(), 'references', 'tier2', sprintf('bfr_%s.mat', slug));
            testCase.assumeTrue(isfile(refFile), ...
                sprintf(['No Tier-2 reference saved yet for BFR method "%s". ', ...
                         'Run regenerate_reference(''confirm'',true,''tier'',''tier2'') to create it after a manual review.'], bfrMethod));

            ref = load(refFile);
            tol = sepiatest.tolerance_for_method('bfr', bfrMethod);
            sepiatest.compare_to_reference(testCase, actual, ref.stats, tol, sprintf('Tier2:BFR:%s', bfrMethod));
        end
    end

end

function methods = get_bfr_methods()
thisFile   = mfilename('fullpath');
testRoot   = fileparts(fileparts(thisFile));
SEPIA_HOME = fileparts(testRoot);
if exist('sepia_universal_variables','file') ~= 2
    addpath(SEPIA_HOME);
end
sepia_universal_variables;
methods = methodBFRName(:)';
end
