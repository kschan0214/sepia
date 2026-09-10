function explore_twopass_masking(varargin)
%% explore_twopass_masking('outputRoot', '/path/to/scratch')
%
% Description: exploratory (not a matlab.unittest / not run in CI) script
% that validates SEPIA's two-pass masking feature on real in-vivo data, and
% sweeps the Magnitude Gradient Field (MGF) two-pass strategy's lambda
% threshold to find a sensible in-vivo default.
%
% Runs across ALL 8 vendor/sequence combinations of the QSM Consensus Paper
% example dataset (https://doi.org/10.1002/mrm.29048), with a fixed
% processing pipeline (per instruction):
%   - Echo phase combination: ROMEO total field calculation, with phase
%     offset correction ON (Wrapper_TotalField_ROMEO.m) - replaces spatial
%     unwrap + temporal echo combination with ROMEO's own 3D+time unwrapping.
%   - Background field removal: VSHARP, default settings.
%   - QSM dipole inversion: "MRI Suscep. Calc." addon, Iterative Tikhonov
%     solver, default settings (alpha=0.05, CG stopping threshold=0.03).
% so two-pass masking is the only isolated variable.
%
% NOTE on inversion method choice: two earlier candidates were tried and
% rejected because they cannot show a two-pass masking effect *by
% construction*, not because of a script bug:
%   - Truncated K-space Division (both SEPIA's built-in TKD and the "MRI
%     Suscep. Calc." addon's TKD solver): TKD.m computes the whole-volume
%     susceptibility map by direct k-space division BEFORE masking - the
%     mask is only a final multiplication (SusceptibilityMap .* Mask),
%     never entering the deconvolution. Confirmed empirically (boundary and
%     core stats were exactly equal between pass-1 and combined at every
%     tested lambda).
%   - Direct Tikhonov (same addon): dirTik.m is the same closed-form
%     k-space-division structure (Kernel = dipole./(dipole.^2+alpha)) with
%     the same post-hoc masking - same null result expected.
% Iterative Tikhonov (iterTik.m) is different: it solves via conjugate
% gradient with the mask baked into the system matrix
% (A = D*W^2*D + alpha*Mask^2), so a refined mask changes the whole-volume
% solution, not just which voxels get zeroed at the end - only this class
% of method can actually exercise two-pass masking's intended benefit.
%
% For each dataset this script:
%   1. Runs the baseline single-pass QSM (isTwoPass='None').
%   2. Sanity-checks all three two-pass refinement strategies once each at
%      their current defaults (Monoexponential decay model, Magnitude
%      Gradient Field, Noise map).
%   3. Sweeps the MGF lambda over a grid and computes, per lambda:
%        - mask retention (Dice / volume ratio vs the pass-1 mask)
%        - an edge-voxel ("boundary shell") streaking-outlier metric
%        - a deep-brain ("core") stability metric
%   4. Saves a .mat file with everything, a markdown summary table, and (for
%      one representative dataset, plus any dataset that stands out as an
%      outlier) slice montage PNGs comparing pass-1 vs combined chi maps.
%
% This script does not change any SEPIA default - it only produces evidence
% to inform that decision. See sepia.documentation's
% docs/method/qsm/Two-pass-masking.rst for the algorithm background.
%
% Usage
% --------------
%   explore_twopass_masking()                          % writes under test/tools/twopass_exploration
%   explore_twopass_masking('outputRoot', '/scratch/x') % writes elsewhere
%   explore_twopass_masking('lambdaGrid', [0.3 0.7 1.5])% override the sweep grid
%
% Kwok-shing Chan
% Date created: 2026-09-09
%
function_name = mfilename;
disp(['--- ' function_name ' ---']);

p = inputParser;
addParameter(p, 'outputRoot', '', @ischar);
addParameter(p, 'lambdaGrid', [0.3 0.5 0.7 1.0 1.5 2.0 3.0], @isnumeric);
addParameter(p, 'datasetFilter', {}, @iscell); % e.g. {'SIEMENS_Monopolar'} to restrict to a subset, for a quick smoke test
parse(p, varargin{:});
lambdaGrid = p.Results.lambdaGrid;
datasetFilter = p.Results.datasetFilter;

%% Setup path (mirrors test/tier3_matrix/TestQSMMatrix.m)
thisFile   = mfilename('fullpath');           % .../test/tools/explore_twopass_masking
testRoot   = fileparts(fileparts(thisFile));  % .../test
SEPIA_HOME = fileparts(testRoot);             % repo root

addpath(SEPIA_HOME);
sepia_addpath; % auto-creates SpecifyToolboxesDirectory.m if missing
addpath(testRoot);

if isempty(p.Results.outputRoot)
    outputRoot = fullfile(testRoot, 'tools', 'twopass_exploration');
else
    outputRoot = p.Results.outputRoot;
end
if ~exist(outputRoot, 'dir'); mkdir(outputRoot); end

%% Check the toolboxes this pipeline needs are configured
run(fullfile(SEPIA_HOME, 'SpecifyToolboxesDirectory.m'));
if isempty(MRITOOLS_HOME) || exist(MRITOOLS_HOME, 'dir') ~= 7
    error('explore_twopass_masking:mritoolsMissing', ...
        ['MRITOOLS_HOME is not configured (or does not exist) in SpecifyToolboxesDirectory.m. ', ...
         'This script uses ROMEO (echoCombMethod) for total field calculation, which needs mritools.']);
end
if isempty(MRISC_HOME) || exist(MRISC_HOME, 'dir') ~= 7
    error('explore_twopass_masking:mriscMissing', ...
        ['MRISC_HOME is not configured (or does not exist) in SpecifyToolboxesDirectory.m. ', ...
         'This script uses the "MRI Suscep. Calc." addon''s Iterative Tikhonov solver as the ', ...
         'fixed dipole-inversion read-out.']);
end

%% Define the 8 QSM Consensus Paper vendor/sequence datasets
% inputDir points at the pre-converted, SEPIA-ready 'derivatives/SEPIA/...'
% directory (single combined multi-echo phase/mag NIfTI + sepia_header.mat)
% - NOT the raw BIDS 'converted/...' source directory. This matches the
% consensus paper's own reference configs
% (QSM_Consensus_Paper_Example_Code/SEPIA_Pipeline_FANSI/SEPIA_*_config.m)
% and also resolves GE's data being shipped as real/imaginary rather than
% phase directly in the raw BIDS directory (the derivatives dir already has
% phase computed). No external mask file is used (isBET=1 instead, matching
% the reference configs - mask_filename='' there). isInvert and
% isEddyCorrect vary per dataset, both taken directly from the matching
% reference config file.
derivRoot = '/autofs/space/symphony_002/users/kwokshing/external_data/QSM_Consensus_Paper_Example_Data_Result_Code_v0.2.1/derivatives/SEPIA';

datasets = struct('label', {}, 'inputDir', {}, 'maskFile', {}, 'isInvert', {}, 'isEddyCorrect', {});
% {vendor, seqDir, isInvert, isEddyCorrect} - isInvert=1 only for GE (both
% readouts); isEddyCorrect=1 only for Bipolar readouts - per
% SEPIA_<VENDOR>_<SEQ>_config.m
combos = { ...
    'SIEMENS', 'Monopolar',         0, 0; ...
    'SIEMENS', 'Bipolar',           0, 1; ...
    'GE',      'Monopolar',         1, 0; ...
    'GE',      'Bipolar',           1, 1; ...
    'PHILIPS', 'Monopolar_CLEAR',   0, 0; ...
    'PHILIPS', 'Bipolar_CLEAR',     0, 1; ...
    'PHILIPS', 'Monopolar_SYNERGY', 0, 0; ...
    'PHILIPS', 'Bipolar_SYNERGY',   0, 1; ...
};
for k = 1:size(combos,1)
    vendor = combos{k,1}; seqDir = combos{k,2};
    ds.label         = sprintf('%s_%s', vendor, seqDir);
    ds.inputDir      = fullfile(derivRoot, vendor, seqDir, 'GRE');
    ds.maskFile      = ''; % isBET=1 - no external mask, matching the reference configs
    ds.isInvert      = combos{k,3};
    ds.isEddyCorrect = combos{k,4};
    datasets(end+1) = ds; %#ok<AGROW>
end

sepia_universal_variables;

% Resumable: if a previous (possibly partial, e.g. crashed mid-way) run
% already saved results in this outputRoot, start from those and only
% (re)process datasets not already present - so a crash on one dataset
% doesn't throw away completed work on the others.
resultsMatFile = fullfile(outputRoot, 'twopass_exploration_results.mat');
if isfile(resultsMatFile)
    prev = load(resultsMatFile, 'results');
    results = prev.results;
    fprintf('Resuming: found %d already-completed dataset(s) in %s\n', numel(results), resultsMatFile);
else
    results = struct('label', {}, 'lambdaGrid', {}, 'strategySanity', {}, 'sweep', {});
end

for iDs = 1:numel(datasets)
    ds = datasets(iDs);

    if ~isempty(datasetFilter) && ~any(strcmp(ds.label, datasetFilter))
        continue
    end
    if any(strcmp(ds.label, {results.label}))
        fprintf('\n=== Dataset %d/%d: %s (already completed, skipping) ===\n', iDs, numel(datasets), ds.label);
        continue
    end

    fprintf('\n=== Dataset %d/%d: %s ===\n', iDs, numel(datasets), ds.label);

    if exist(ds.inputDir, 'dir') ~= 7
        warning('explore_twopass_masking:datasetMissing', ...
            'Dataset "%s" not found on this machine (inputDir missing) - skipping.', ds.label);
        continue
    end

    dsOutRoot = fullfile(outputRoot, ds.label);
    if ~exist(dsOutRoot, 'dir'); mkdir(dsOutRoot); end

    try

    % Matches the consensus-paper's own reference config exactly
    % (QSM_Consensus_Paper_Example_Code/SEPIA_Pipeline_FANSI/
    %  SEPIA_<VENDOR>_<SEQ>_config.m), except the QSM inversion method
    % (that reference uses FANSI; here it's fixed to MRI Suscep. Calc. /
    % Iterative Tikhonov to isolate the two-pass masking effect - see file
    % header for why). isInvert and isEddyCorrect come from ds (per-dataset,
    % taken from the matching reference config file).
    baseAlgorParam = struct();
    baseAlgorParam.general.isBET               = 1; % no external mask - BET run internally
    baseAlgorParam.general.fractional_threshold = 0.5;
    baseAlgorParam.general.gradient_threshold   = 0;
    baseAlgorParam.general.isInvert             = ds.isInvert;

    % Echo phase combination: ROMEO total field calculation, phase offset
    % correction ON. ROMEO computes totalField directly from the wrapped
    % multi-echo phase (its own 3D+time unwrapping), so unwrap.unwrapMethod
    % (spatial-unwrap-per-echo) is not used in this path.
    baseAlgorParam.unwrap.echoCombMethod       = 'ROMEO total field calculation';
    baseAlgorParam.unwrap.offsetCorrect        = 'On';
    baseAlgorParam.unwrap.mask                 = 'SEPIA mask';
    baseAlgorParam.unwrap.qualitymaskThreshold = 0.5;
    baseAlgorParam.unwrap.useRomeoMask         = false;
    baseAlgorParam.unwrap.isEddyCorrect        = ds.isEddyCorrect;
    baseAlgorParam.unwrap.isSaveUnwrappedEcho  = 0;
    baseAlgorParam.unwrap.excludeMaskThreshold = 0.3;
    baseAlgorParam.unwrap.excludeMethod        = 'Weighting map';

    % Background field removal: VSHARP, matching the consensus-paper
    % example config's radius range (12mm down to 1mm), no polynomial
    % refinement step.
    baseAlgorParam.bfr.method             = 'VSHARP';
    baseAlgorParam.bfr.refine_method      = 'None';
    baseAlgorParam.bfr.refine_order       = 2;
    baseAlgorParam.bfr.erode_before_radius = 1;
    baseAlgorParam.bfr.erode_radius       = 0;
    baseAlgorParam.bfr.radius             = 12:-1:1;

    % QSM dipole inversion: "MRI Suscep. Calc." addon, Iterative Tikhonov
    % solver, default settings (alpha=0.05, CG stopping threshold=0.03) -
    % see file header for why this solver (and not TKD/Direct Tikhonov)
    % is required for two-pass masking to have any effect.
    baseAlgorParam.qsm.method = 'MRI Suscep. Calc.';
    baseAlgorParam.qsm.solver = 'Iterative Tikhonov';

    %% 1. Baseline (single pass)
    outPrefix = fullfile(dsOutRoot, 'baseline', 'sepia');
    algorParam = baseAlgorParam;
    algorParam.qsm.isTwoPass = 'None';
    sepiaIO(ds.inputDir, outPrefix, ds.maskFile, algorParam);
    addpath(testRoot); % sepiaIO->sepia_addpath strips test/ off the path

    baselineChi  = double(load_nii_img_only([outPrefix '_Chimap.nii.gz']));
    maskPass1    = double(load_nii_img_only([outPrefix '_mask_QSM.nii.gz'])) > 0;

    %% 2. Strategy sanity-check (all three, at their current defaults)
    strategySanity = struct('strategy', {}, 'volumeRatio', {}, 'dice', {}, 'ranOk', {});
    for iStrat = 2:4 % methodTwoPassName{1}='None'; {2,3,4} are the three strategies
        strategy = methodTwoPassName{iStrat};
        row.strategy = strategy;
        row.ranOk = true;
        try
            stratPrefix = fullfile(dsOutRoot, ['strategy_' matlab.lang.makeValidName(strategy)], 'sepia');
            algorParam = baseAlgorParam;
            algorParam.qsm.isTwoPass      = strategy;
            algorParam.qsm.twopass_lambda = 0.7; % current global default
            sepiaIO(ds.inputDir, stratPrefix, ds.maskFile, algorParam);
            addpath(testRoot);

            maskRefined = double(load_nii_img_only([stratPrefix '_mask_QSM-2pass.nii.gz'])) > 0;
            row.volumeRatio = nnz(maskRefined) / max(nnz(maskPass1), 1);
            row.dice        = dice_coeff(maskPass1, maskRefined);
        catch ME
            warning('explore_twopass_masking:strategyFailed', ...
                'Strategy "%s" failed on dataset "%s": %s', strategy, ds.label, ME.message);
            row.ranOk = false;
            row.volumeRatio = NaN;
            row.dice = NaN;
        end
        strategySanity(end+1) = row; %#ok<AGROW>
    end

    %% 3. MGF lambda sweep
    % "boundary" is redefined per-lambda below as the band of voxels RETAINED
    % in the refined mask but adjacent to voxels the refinement excluded -
    % this is where a two-pass benefit (streaking reduction) can actually
    % show up. (A fixed shell at the *original* mask edge would trivially
    % show zero difference, since excluded edge voxels are never touched by
    % pass 2 - chi_combined keeps the pass-1 value there by construction.)
    deepCore = imerode(maskPass1, strel('sphere', 9));

    sweep = struct('lambda', {}, 'volumeRatio', {}, 'dice', {}, 'nBoundaryVoxels', {}, ...
        'boundaryPass1', {}, 'boundaryCombined', {}, 'corePass1', {}, 'coreCombined', {});
    for iL = 1:numel(lambdaGrid)
        lambda = lambdaGrid(iL);
        fprintf('  MGF lambda = %.2f ...\n', lambda);
        lamPrefix = fullfile(dsOutRoot, sprintf('mgf_lambda_%03d', round(lambda*100)), 'sepia');

        algorParam = baseAlgorParam;
        algorParam.qsm.isTwoPass      = methodTwoPassName{3}; % Magnitude Gradient Field
        algorParam.qsm.twopass_lambda = lambda;
        sepiaIO(ds.inputDir, lamPrefix, ds.maskFile, algorParam);
        addpath(testRoot);

        maskRefined  = double(load_nii_img_only([lamPrefix '_mask_QSM-2pass.nii.gz'])) > 0;
        chiCombined  = double(load_nii_img_only([lamPrefix '_Chimap.nii.gz']));

        excludedVoxels = maskPass1 & ~maskRefined;
        nearExcludedBand = imdilate(excludedVoxels, strel('sphere', 3)) & maskRefined;

        clear srow
        srow.lambda      = lambda;
        srow.volumeRatio = nnz(maskRefined) / max(nnz(maskPass1), 1);
        srow.dice        = dice_coeff(maskPass1, maskRefined);
        srow.nBoundaryVoxels = nnz(nearExcludedBand);

        srow.boundaryPass1    = sepiatest.mask_stats(baselineChi, nearExcludedBand);
        srow.boundaryCombined = sepiatest.mask_stats(chiCombined, nearExcludedBand);
        srow.corePass1        = sepiatest.mask_stats(baselineChi, deepCore);
        srow.coreCombined     = sepiatest.mask_stats(chiCombined, deepCore);

        sweep(end+1) = srow; %#ok<AGROW>

        % visual QC for a representative dataset only, at a subset of lambdas
        if strcmp(ds.label, 'SIEMENS_Monopolar') && any(abs(lambda - [0.3 0.7 1.5 3.0]) < 1e-6)
            save_slice_montage(baselineChi, chiCombined, maskRefined, ...
                fullfile(outputRoot, sprintf('montage_%s_lambda%.1f.png', ds.label, lambda)), ...
                sprintf('%s, MGF lambda=%.1f', strrep(ds.label,'_',' '), lambda));
        end
    end

    results(end+1).label          = ds.label; %#ok<AGROW>
    results(end).lambdaGrid       = lambdaGrid;
    results(end).strategySanity   = strategySanity;
    results(end).sweep            = sweep;

    % save incrementally so a crash on a LATER dataset never loses this one
    save(resultsMatFile, 'results', 'lambdaGrid', '-v7.3');
    write_markdown_report(results, fullfile(outputRoot, 'twopass_exploration_report.md'));

    catch ME
        warning('explore_twopass_masking:datasetFailed', ...
            'Dataset "%s" failed and was skipped: %s', ds.label, ME.message);
        addpath(testRoot); % defensive: an error inside sepiaIO may have left test/ stripped off the path
    end
end

%% Save everything + write markdown summary
save(fullfile(outputRoot, 'twopass_exploration_results.mat'), 'results', 'lambdaGrid', '-v7.3');
write_markdown_report(results, fullfile(outputRoot, 'twopass_exploration_report.md'));

fprintf('\nDone. Results saved under: %s\n', outputRoot);

end

%% ------------------------------------------------------------------
function d = dice_coeff(maskA, maskB)
maskA = maskA > 0; maskB = maskB > 0;
denom = nnz(maskA) + nnz(maskB);
if denom == 0
    d = NaN;
else
    d = 2 * nnz(maskA & maskB) / denom;
end
end

%% ------------------------------------------------------------------
function save_slice_montage(chiPass1, chiCombined, maskRefined, outFile, titleStr)
sz = size(chiPass1);
sliceIdx = round(linspace(round(sz(3)*0.3), round(sz(3)*0.7), 4));

fig = figure('Visible', 'off', 'Position', [0 0 1200 700]);
tiledlayout(2, numel(sliceIdx), 'TileSpacing', 'compact', 'Padding', 'compact');
climVal = [-0.15 0.15];
for i = 1:numel(sliceIdx)
    nexttile;
    imagesc(rot90(chiPass1(:,:,sliceIdx(i))), climVal); axis image off; colormap(gca, gray);
    if i == 1; ylabel('pass-1 (single pass)'); end
    title(sprintf('slice %d', sliceIdx(i)));
end
for i = 1:numel(sliceIdx)
    nexttile;
    imagesc(rot90(chiCombined(:,:,sliceIdx(i))), climVal); axis image off; colormap(gca, gray);
    if i == 1; ylabel('two-pass combined'); end
end
sgtitle(titleStr, 'Interpreter', 'none');
exportgraphics(fig, outFile, 'Resolution', 120);
close(fig);
end

%% ------------------------------------------------------------------
function write_markdown_report(results, outFile)
fid = fopen(outFile, 'w');
fprintf(fid, '# Two-pass masking exploration report\n\n');
fprintf(fid, ['Fixed pipeline: ROMEO total field calculation (offset correction on) + VSHARP ', ...
              '(radius 12:-1:1, no polynomial refine) + MRI Suscep. Calc. Iterative Tikhonov ', ...
              '(alpha=0.05, tol=0.03), isBET=1, matching the consensus-paper reference configs.\n\n']);

for iDs = 1:numel(results)
    r = results(iDs);
    fprintf(fid, '## %s\n\n', strrep(r.label, '_', ' '));

    fprintf(fid, '### Strategy sanity check (lambda=0.7 / default threshold)\n\n');
    fprintf(fid, '| strategy | volume ratio | dice | ran ok |\n|---|---|---|---|\n');
    for s = r.strategySanity
        fprintf(fid, '| %s | %.3f | %.3f | %d |\n', s.strategy, s.volumeRatio, s.dice, s.ranOk);
    end
    fprintf(fid, '\n');

    fprintf(fid, '### MGF lambda sweep\n\n');
    fprintf(fid, ['"boundary" = voxels retained in the refined mask but adjacent (within 3 vox) ', ...
                  'to a voxel the refinement excluded - where a two-pass benefit should show up. ', ...
                  '"core" = mask eroded by 9 vox, a deep-brain stability check (should stay ~unchanged).\n\n']);
    fprintf(fid, ['| lambda | mask vol ratio | dice | n boundary vox | boundary std (pass1->combined) | ', ...
                  'boundary p1/p99 (pass1->combined) | core mean (pass1->combined) |\n|---|---|---|---|---|---|---|\n']);
    for s = r.sweep
        fprintf(fid, '| %.2f | %.3f | %.3f | %d | %.4f -> %.4f | [%.4f,%.4f] -> [%.4f,%.4f] | %.4f -> %.4f |\n', ...
            s.lambda, s.volumeRatio, s.dice, s.nBoundaryVoxels, ...
            s.boundaryPass1.std, s.boundaryCombined.std, ...
            s.boundaryPass1.p1, s.boundaryPass1.p99, s.boundaryCombined.p1, s.boundaryCombined.p99, ...
            s.corePass1.mean, s.coreCombined.mean);
    end
    fprintf(fid, '\n');
end

fclose(fid);
end
