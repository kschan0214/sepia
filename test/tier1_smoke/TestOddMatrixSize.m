%% TestOddMatrixSize - Tier 1 toolbox-free odd-matrix-size test
%
% SEPIA pads odd-sized matrix dimensions to even before FFT/k-space-based
% processing (estimateTotalField, BackgroundRemovalMacro, QSMMacro - see
% utils/zeropad_odd_dimension.m) and crops back to the original size
% afterward. This test runs the full one-stop pipeline on a synthetic
% phantom generated at several odd/even matrix sizes and checks that every
% output comes back at exactly the original (unpadded) size with no
% NaN/Inf - i.e. that the pad/crop round trip is correct.
%
% It also exercises get_variable_from_headerAndExtraData.m's disk-loaded
% auxiliary-data padding path via Wrapper_QSM_iLSQR.m's 'weights' input,
% which the default TestSmokePhantom pipeline never supplies at all.
%
% Only covers what's reachable without an external toolbox: estimateTotalField
% and BackgroundRemovalMacro/QSMMacro pad unconditionally regardless of
% method, so unwrap='None' + bfr='VSHARP' + qsm='TKD'/'iLSQR' already
% exercises them. UnwrapPhaseMacro's OWN pad/crop only triggers for
% Laplacian-based unwrap methods, which require MEDI/STI Suite - not
% reachable here (see test/README.md and the plan this test came from).
%
classdef TestOddMatrixSize < matlab.unittest.TestCase

    properties (TestParameter)
        matrixSize = {[32 32 24], [31 32 24], [32 31 24], [32 32 25], [31 31 25]};
    end

    methods (TestClassSetup)
        function addSepiaToPath(testCase)
            % see TestSmokePhantom.m for why test/ must be re-added after
            % sepia_addpath (which does rmpath(genpath(SEPIA_HOME)))
            thisFile   = mfilename('fullpath');          % .../test/tier1_smoke/TestOddMatrixSize
            testRoot   = fileparts(fileparts(thisFile)); % .../test
            SEPIA_HOME = fileparts(testRoot);            % repo root

            addpath(SEPIA_HOME);
            sepia_addpath;

            addpath(testRoot);
            addpath(fullfile(testRoot, 'phantom'));
        end
    end

    methods (Test)
        function testOneStopPipelineShapeRoundTrip(testCase, matrixSize)
            import matlab.unittest.fixtures.TemporaryFolderFixture

            fixture = testCase.applyFixture(TemporaryFolderFixture);
            phantomPaths = generate_synthetic_phantom(fullfile(fixture.Folder, 'phantom'), matrixSize);

            outFixture   = testCase.applyFixture(TemporaryFolderFixture);
            outputPrefix = fullfile(outFixture.Folder, 'sepia');

            algorParam = struct();
            algorParam.general.isInvert = false;
            algorParam.general.isBET    = false;
            algorParam.unwrap.unwrapMethod   = 'None';
            algorParam.unwrap.echoCombMethod = 'Optimum weights';
            algorParam.bfr.method  = 'VSHARP';
            algorParam.qsm.method  = 'TKD';

            sepiaIO(phantomPaths.input, outputPrefix, phantomPaths.mask, algorParam);

            % sepiaIO strips test/, this file's own tier1_smoke/ folder,
            % and test/phantom/ back off the path via sepia_addpath - all
            % three must be re-added, or (respectively) sepiatest.*, this
            % class's own dispatch, and the next parameterized case's call
            % to generate_synthetic_phantom break.
            addpath(fileparts(fileparts(mfilename('fullpath'))));
            addpath(fileparts(mfilename('fullpath')));
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'phantom'));

            sepiatest.assert_output_contract(testCase, outputPrefix, '.nii.gz', ...
                { {'Chimap', true}, {'localfield', true}, {'fieldmap', true} }, ...
                phantomPaths.mask);
        end

        function testAuxiliaryDataLoadingPadsCorrectly(testCase, matrixSize)
            import matlab.unittest.fixtures.TemporaryFolderFixture

            fixture = testCase.applyFixture(TemporaryFolderFixture);
            phantomPaths = generate_synthetic_phantom(fullfile(fixture.Folder, 'phantom'), matrixSize);

            outFixture   = testCase.applyFixture(TemporaryFolderFixture);
            outputPrefix = fullfile(outFixture.Folder, 'sepia');

            % supply a weights file (reusing the phantom's own 3D mask
            % NIfTI as a stand-in - weights must be 3D, and this test is
            % about the load/pad mechanism in
            % get_variable_from_headerAndExtraData.m working correctly at
            % odd matrix sizes, not about weighting-algorithm accuracy)
            input = phantomPaths.input;
            input(3).name = phantomPaths.mask;

            algorParam = struct();
            algorParam.general.isInvert = false;
            algorParam.general.isBET    = false;
            algorParam.unwrap.unwrapMethod   = 'None';
            algorParam.unwrap.echoCombMethod = 'Optimum weights';
            algorParam.bfr.method  = 'VSHARP';
            algorParam.qsm.method  = 'iLSQR';

            sepiaIO(input, outputPrefix, phantomPaths.mask, algorParam);

            addpath(fileparts(fileparts(mfilename('fullpath'))));
            addpath(fileparts(mfilename('fullpath')));
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'phantom'));

            sepiatest.assert_output_contract(testCase, outputPrefix, '.nii.gz', ...
                { {'Chimap', true} }, ...
                phantomPaths.mask);
        end
    end

end
