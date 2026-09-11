%% TestUnwrapMatrixPhantom - Tier 2 synthetic-phantom phase-unwrapping method regression matrix
%
% Mirrors TestQSMMatrixPhantom.m: runs the one-stop SEPIA pipeline on the
% synthetic phantom once per phase-unwrapping method in methodUnwrapName
% (sepia_universal_variables.m) AND per matrix size (the same 5 odd/even
% sizes TestOddMatrixSize.m uses in Tier 1), holding background field
% removal ('VSHARP') and QSM dipole inversion ('TKD', both toolbox-free)
% fixed so the unwrap method is the only isolated variable. '3D best path'
% is deliberately excluded (see sepiatest.unwrap_method_toolbox_key.m).
%
% This is the concrete place UnwrapPhaseMacro.m's own pad/crop round trip
% (Laplacian-based methods only - see utils/zeropad_odd_dimension.m) gets
% exercised at odd matrix sizes: it needs MEDI/STI Suite, so it's
% unreachable in Tier 1's toolbox-free TestOddMatrixSize.m.
%
classdef TestUnwrapMatrixPhantom < matlab.unittest.TestCase

    properties (TestParameter)
        unwrapMethod = get_unwrap_methods();
        matrixSize   = {[32 32 24], [31 32 24], [32 31 24], [32 32 25], [31 31 25]};
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
        function testUnwrapMethodMatchesReference(testCase, unwrapMethod, matrixSize)
            import matlab.unittest.fixtures.TemporaryFolderFixture

            addpath(fileparts(fileparts(mfilename('fullpath'))));
            toolboxKey = sepiatest.unwrap_method_toolbox_key(unwrapMethod);
            sepiatest.assume_toolbox_available(testCase, testCase.Toolboxes, toolboxKey, unwrapMethod);

            fixture = testCase.applyFixture(TemporaryFolderFixture);
            phantomPaths = generate_synthetic_phantom(fullfile(fixture.Folder, 'phantom'), matrixSize);

            outFixture   = testCase.applyFixture(TemporaryFolderFixture);
            outputPrefix = fullfile(outFixture.Folder, 'sepia');

            algorParam = struct();
            algorParam.general.isInvert = false;
            algorParam.general.isBET    = false;
            algorParam.unwrap.unwrapMethod   = unwrapMethod;
            algorParam.unwrap.echoCombMethod = 'Optimum weights';
            algorParam.bfr.method = 'VSHARP';
            algorParam.qsm.method = 'TKD';

            sepiaIO(phantomPaths.input, outputPrefix, phantomPaths.mask, algorParam);

            addpath(fileparts(fileparts(mfilename('fullpath'))));
            addpath(fileparts(mfilename('fullpath')));
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'phantom'));

            suffix = '.nii.gz';
            slug   = matlab.lang.makeValidName(unwrapMethod);

            sepiatest.assert_output_contract(testCase, outputPrefix, suffix, ...
                { {'fieldmap', true}, {'Chimap', true} }, phantomPaths.mask);

            if ~isequal(matrixSize, [32 32 24])
                return
            end

            maskImg   = load_nii_img_only(phantomPaths.mask) > 0;
            fieldImg  = load_nii_img_only([outputPrefix '_fieldmap' suffix]);
            actual    = sepiatest.mask_stats(fieldImg, maskImg);

            refFile = fullfile(sepiatest.test_root(), 'references', 'tier2', sprintf('unwrap_%s.mat', slug));
            testCase.assumeTrue(isfile(refFile), ...
                sprintf(['No Tier-2 reference saved yet for unwrap method "%s". ', ...
                         'Run regenerate_reference(''confirm'',true,''tier'',''tier2'') to create it after a manual review.'], unwrapMethod));

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
