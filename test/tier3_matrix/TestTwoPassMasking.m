%% TestTwoPassMasking - Tier 3 real-dataset two-pass masking regression matrix
%
% Runs the one-stop SEPIA pipeline on a real dataset (path supplied via
% sepiatest.get_real_dataset(), never hardcoded) once per two-pass mask
% refinement strategy in methodTwoPassName (sepia_universal_variables.m,
% excluding 'None'), holding phase unwrapping ('None'), background field
% removal ('VSHARP', the built-in variant) and the QSM dipole-inversion
% method (FANSI) fixed so the refinement strategy is the only isolated
% variable.
%
% FANSI is used as the fixed QSM method (not TKD, unlike the other Tier 3
% matrices) because two-pass masking has NO effect on TKD/Direct
% Tikhonov's output by construction - both apply the mask only as a final
% multiplication after a closed-form/direct k-space inversion, so the mask
% never enters the deconvolution itself (see
% sepia.documentation/docs/method/qsm/Two-pass-masking.rst and
% utils/qsm_method_uses_mask_in_inversion.m). FANSI's mask enters its
% data-fidelity term directly, so a refined mask can actually change the
% reconstructed values - this is the class of method the exploratory sweep
% in test/tools/explore_twopass_masking.m used to validate the feature.
%
% Each row is skipped (not failed) when:
%   - no real dataset is configured on this machine, or
%   - FANSI isn't installed, or
%   - no reference has been saved yet for that strategy
%
classdef TestTwoPassMasking < matlab.unittest.TestCase

    properties (TestParameter)
        twoPassStrategy = get_two_pass_strategies();
    end

    properties
        Dataset
        Toolboxes
    end

    methods (TestClassSetup)
        function setupPathAndData(testCase)
            thisFile   = mfilename('fullpath');          % .../test/tier3_matrix/TestTwoPassMasking
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
        function testTwoPassMatchesReference(testCase, twoPassStrategy)
            import matlab.unittest.fixtures.TemporaryFolderFixture

            addpath(fileparts(fileparts(mfilename('fullpath')))); % defensive: ensure test/ is still on path
            sepiatest.assume_toolbox_available(testCase, testCase.Toolboxes, 'FANSI', 'FANSI');

            fixture = testCase.applyFixture(TemporaryFolderFixture);
            outputPrefix = fullfile(fixture.Folder, 'sepia');

            algorParam = struct();
            algorParam.general.isInvert = false;
            algorParam.general.isBET    = false;
            algorParam.unwrap.unwrapMethod   = 'None';
            algorParam.unwrap.echoCombMethod = 'Optimum weights';
            algorParam.bfr.method = 'VSHARP';
            algorParam.qsm.method = 'FANSI';
            algorParam.qsm.isGPU  = false; % CPU-only for determinism
            algorParam.qsm.isTwoPass = twoPassStrategy;

            sepiaIO(testCase.Dataset.inputDir, outputPrefix, testCase.Dataset.maskFile, algorParam);

            % sepiaIO/SepiaIOWrapper calls sepia_addpath internally, which
            % strips test/ back off the path - re-add before using
            % sepiatest.* below.
            addpath(fileparts(fileparts(mfilename('fullpath'))));

            suffix = '.nii.gz';
            slug   = matlab.lang.makeValidName(twoPassStrategy);

            sepiatest.assert_output_contract(testCase, outputPrefix, suffix, ...
                { {'Chimap', true}, {'desc-firstpass_Chimap', true}, {'desc-secondpass_Chimap', true}, ...
                  {'mask_QSM-2pass', true} }, ...
                testCase.Dataset.maskFile);

            maskImg = load_nii_img_only(testCase.Dataset.maskFile) > 0;

            actual = struct();
            actual.combined = sepiatest.mask_stats(load_nii_img_only([outputPrefix '_Chimap' suffix]), maskImg);
            actual.pass1    = sepiatest.mask_stats(load_nii_img_only([outputPrefix '_desc-firstpass_Chimap' suffix]), maskImg);
            actual.pass2    = sepiatest.mask_stats(load_nii_img_only([outputPrefix '_desc-secondpass_Chimap' suffix]), maskImg);
            refinedMask     = load_nii_img_only([outputPrefix '_mask_QSM-2pass' suffix]) > 0;
            actual.maskVolumeRatio = nnz(refinedMask) / max(nnz(maskImg),1);

            refFile = fullfile(sepiatest.test_root(), 'references', 'tier3', sprintf('twopass_%s.mat', slug));
            testCase.assumeTrue(isfile(refFile), ...
                sprintf(['No Tier-3 reference saved yet for two-pass strategy "%s". ', ...
                         'Run regenerate_reference(''confirm'',true,''tier'',''tier3'') to create it after a manual review.'], twoPassStrategy));

            ref = load(refFile);
            tol = sepiatest.tolerance_for_method('twopass', twoPassStrategy);
            label = sprintf('Tier3:TwoPass:%s', twoPassStrategy);
            sepiatest.compare_to_reference(testCase, actual.combined, ref.stats.combined, tol, [label ':combined']);
            sepiatest.compare_to_reference(testCase, actual.pass1,    ref.stats.pass1,    tol, [label ':pass1']);
            sepiatest.compare_to_reference(testCase, actual.pass2,    ref.stats.pass2,    tol, [label ':pass2']);
            testCase.verifyEqual(actual.maskVolumeRatio, ref.stats.maskVolumeRatio, ...
                'RelTol', tol.relTol, 'AbsTol', tol.absTol, ...
                sprintf('%s: refined mask volume ratio outside tolerance (actual=%.4f, reference=%.4f)', ...
                        label, actual.maskVolumeRatio, ref.stats.maskVolumeRatio));
        end
    end

end

%% pull the strategy list from SEPIA itself rather than hardcoding it here
function strategies = get_two_pass_strategies()

thisFile   = mfilename('fullpath');           % .../test/tier3_matrix/TestTwoPassMasking
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
strategies = methodTwoPassName(2:end); % exclude 'None' - a single-pass baseline isn't a two-pass row
strategies = strategies(:)';

end
