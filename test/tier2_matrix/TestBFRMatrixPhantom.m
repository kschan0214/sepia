%% TestBFRMatrixPhantom - Tier 2 synthetic-phantom background-field-removal method regression matrix
%
% Mirrors TestQSMMatrixPhantom.m: runs the one-stop SEPIA pipeline on the
% synthetic phantom once per background-field-removal method in
% methodBFRName (sepia_universal_variables.m) AND per matrix size (the
% same 5 odd/even sizes TestOddMatrixSize.m uses in Tier 1), holding phase
% unwrapping ('None') and QSM dipole inversion ('TKD', toolbox-free) fixed
% so the BFR method is the only isolated variable. Only the baseline
% matrix size ([32 32 24]) is compared numerically against a saved
% reference; the odd-size cases only assert output shape/no-NaN.
%
classdef TestBFRMatrixPhantom < matlab.unittest.TestCase

    properties (TestParameter)
        bfrMethod  = get_bfr_methods();
        matrixSize = {[32 32 24], [31 32 24], [32 31 24], [32 32 25], [31 31 25]};
    end

    properties
        Toolboxes
    end

    methods (TestClassSetup)
        function setupPathAndToolboxes(testCase)
            thisFile   = mfilename('fullpath');
            testRoot   = fileparts(fileparts(thisFile));
            SEPIA_HOME = fileparts(testRoot);

            addpath(SEPIA_HOME);
            sepia_addpath;
            addpath(testRoot);
            addpath(fullfile(testRoot, 'phantom'));

            testCase.Toolboxes = sepiatest.discover_toolboxes();
            addpath(testRoot); addpath(fullfile(testRoot, 'phantom'));
        end
    end

    methods (Test)
        function testBFRMethodMatchesReference(testCase, bfrMethod, matrixSize)
            import matlab.unittest.fixtures.TemporaryFolderFixture

            addpath(fileparts(fileparts(mfilename('fullpath'))));
            toolboxKey = sepiatest.bfr_method_toolbox_key(bfrMethod);
            sepiatest.assume_toolbox_available(testCase, testCase.Toolboxes, toolboxKey, bfrMethod);

            fixture = testCase.applyFixture(TemporaryFolderFixture);
            phantomPaths = generate_synthetic_phantom(fullfile(fixture.Folder, 'phantom'), matrixSize);

            outFixture   = testCase.applyFixture(TemporaryFolderFixture);
            outputPrefix = fullfile(outFixture.Folder, 'sepia');

            algorParam = struct();
            algorParam.general.isInvert = false;
            algorParam.general.isBET    = false;
            algorParam.unwrap.unwrapMethod   = 'None';
            algorParam.unwrap.echoCombMethod = 'Optimum weights';
            algorParam.bfr.method = bfrMethod;
            algorParam.qsm.method = 'TKD';

            sepiaIO(phantomPaths.input, outputPrefix, phantomPaths.mask, algorParam);

            addpath(fileparts(fileparts(mfilename('fullpath'))));
            addpath(fileparts(mfilename('fullpath')));
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'phantom'));

            suffix = '.nii.gz';
            slug   = matlab.lang.makeValidName(bfrMethod);

            sepiatest.assert_output_contract(testCase, outputPrefix, suffix, ...
                { {'localfield', true}, {'Chimap', true} }, phantomPaths.mask);

            if ~isequal(matrixSize, [32 32 24])
                return
            end

            maskImg = load_nii_img_only(phantomPaths.mask) > 0;
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
methods = methodBFRName(:)';
end
