%% outputFileList = construct_output_filename(outputDir, ouputPrefix, outputSufffix)
%
% Input
% --------------
% outputDir     : output directory
% ouputPrefix   : output filename prefix
%
% Output
% --------------
% outputFileList: structure contains all output filenames
%
% Kwok-shing Chan @ DCCN
% kwokshing.chan@donders.ru.nl
% Date created: 11 August 2021 (v1.0)
% Date modified:
%
%
function [outputFileList,ouputPrefix] = construct_output_filename(outputDir, ouputPrefix, algorParam, outputSuffix)

% This function assigns its own 'desc-' label to (some of) its outputs to
% distinguish output types (e.g. desc-rad, desc-unwrapped,
% desc-paramagnetic). BIDS filenames can only contain one 'desc-' entity,
% so if ouputPrefix already carries one, merge it with SEPIA's own
% label(s) into a single camelCase-joined value, in the order processing
% actually happens, instead of ending up with two 'desc-' entities.
[ouputPrefix, existingDesc] = extract_and_remove_desc_entity(ouputPrefix);

% phase related: independent, one-off outputs, each merged only with
% whatever desc- was already in the prefix (not with one another)
outputFileList.phaseRadian      = fullfile(outputDir, [ouputPrefix 'part-phase_desc-' join_desc_label(existingDesc,'rad')          outputSuffix]);
outputFileList.phaseReversed    = fullfile(outputDir, [ouputPrefix 'part-phase_desc-' join_desc_label(existingDesc,'reverse')      outputSuffix]);
outputFileList.phaseEddyCorr    = fullfile(outputDir, [ouputPrefix 'part-phase_desc-' join_desc_label(existingDesc,'bipolarcorr')  outputSuffix]);
outputFileList.unwrappedPhase   = fullfile(outputDir, [ouputPrefix 'part-phase_desc-' join_desc_label(existingDesc,'unwrapped')    outputSuffix]);

% standard output
outputFileList.totalField       = fullfile(outputDir, [ouputPrefix 'fieldmap' outputSuffix]);
outputFileList.localField       = fullfile(outputDir, [ouputPrefix 'localfield' outputSuffix]);
outputFileList.QSM              = fullfile(outputDir, [ouputPrefix 'Chimap' outputSuffix]);
outputFileList.QSMpara          = fullfile(outputDir, [ouputPrefix 'desc-' join_desc_label(existingDesc,'paramagnetic') '_Chimap' outputSuffix]);
outputFileList.QSMdia           = fullfile(outputDir, [ouputPrefix 'desc-' join_desc_label(existingDesc,'diamagnetic') '_Chimap' outputSuffix]);

% use for regularisation
outputFileList.weights          = fullfile(outputDir, [ouputPrefix 'weights' outputSuffix]);
outputFileList.fieldmapSD       = fullfile(outputDir, [ouputPrefix 'noisesd' outputSuffix]);
outputFileList.relativeResidual	= fullfile(outputDir, [ouputPrefix 'relativeresidual' outputSuffix]);
outputFileList.relativeResidualWeights	= fullfile(outputDir, [ouputPrefix 'relativeresidualweights' outputSuffix]);

% derived masks
outputFileList.maskBrain        = fullfile(outputDir, [ouputPrefix 'mask_brain' outputSuffix]);
outputFileList.maskReliable     = fullfile(outputDir, [ouputPrefix 'mask_reliable' outputSuffix]);
outputFileList.maskLocalField 	= fullfile(outputDir, [ouputPrefix 'mask_localfield' outputSuffix]);
outputFileList.maskQSM          = fullfile(outputDir, [ouputPrefix 'mask_QSM' outputSuffix]);
outputFileList.maskQSM2pass     = fullfile(outputDir, [ouputPrefix 'mask_QSM-2pass' outputSuffix]);
outputFileList.maskRef      	= fullfile(outputDir, [ouputPrefix 'mask_referenceregion' outputSuffix]);

% R2*
outputFileList.r2s              = fullfile(outputDir, [ouputPrefix 'R2starmap' outputSuffix]);
outputFileList.t2s              = fullfile(outputDir, [ouputPrefix 'T2starmap' outputSuffix]);
outputFileList.s0               = fullfile(outputDir, [ouputPrefix 'S0map' outputSuffix]);

% misc
outputFileList.phase_bipolar    = fullfile(outputDir, [ouputPrefix 'bipolar_phase' outputSuffix]);
outputFileList.optimalCombinedMagnitude = fullfile(outputDir, [ouputPrefix 'part-mag_desc-' join_desc_label(existingDesc,'optimalcombined') outputSuffix]);

% denoising and upsampling are sequential steps on the SAME data
% (upsampling runs on the already-denoised data when both are enabled),
% so their desc- labels are chained in that processing order rather than
% each being merged with existingDesc independently.
descChain = existingDesc;

if algorParam.general.isDenoise
descChain = join_desc_label(descChain, 'denoised');
outputFileList.magDenoise       = fullfile(outputDir, [ouputPrefix 'part-mag_desc-'   descChain outputSuffix]);
outputFileList.phaseDenoise     = fullfile(outputDir, [ouputPrefix 'part-phase_desc-' descChain outputSuffix]);
outputFileList.sigma            = fullfile(outputDir, [ouputPrefix 'sigma' outputSuffix]);
outputFileList.snrgain          = fullfile(outputDir, [ouputPrefix 'SNRgain' outputSuffix]);
outputFileList.P                = fullfile(outputDir, [ouputPrefix 'P' outputSuffix]);
end

if algorParam.general.isUpsample
descChain = join_desc_label(descChain, 'upsampled');
outputFileList.magUpsample      = fullfile(outputDir, [ouputPrefix 'part-mag_desc-'   descChain outputSuffix]);
outputFileList.phaseUpsample    = fullfile(outputDir, [ouputPrefix 'part-phase_desc-' descChain outputSuffix]);
outputFileList.maskUpsample     = fullfile(outputDir, [ouputPrefix 'mask_upsampled' outputSuffix]);
outputFileList.sepiaHeaderUpsample = fullfile(outputDir, [ouputPrefix 'sepia_header.mat']);
end

end

%% Remove a 'desc-<label>_' entity from a BIDS-style prefix, if present
function [prefix, descLabel] = extract_and_remove_desc_entity(prefix)

tok = regexp(prefix, 'desc-([A-Za-z0-9]+)_', 'tokens', 'once');
if isempty(tok)
    descLabel = '';
else
    descLabel = tok{1};
    prefix    = regexprep(prefix, 'desc-[A-Za-z0-9]+_', '', 'once');
    warning('construct_output_filename:mergedDescEntity', ...
        ['The output prefix already contains a ''desc-%s'' entity. Since BIDS filenames can only ', ...
         'contain one ''desc-'' entity, it will be merged with SEPIA''s own output-type label(s) ', ...
         '(e.g. ''desc-%sRad'').'], descLabel, descLabel);
end

end

%% Combine an existing desc- label (if any) with a new one into a single BIDS-valid, camelCase value
function combined = join_desc_label(existingDesc, newLabel)

if isempty(existingDesc)
    combined = newLabel;
else
    combined = [existingDesc, upper(newLabel(1)), newLabel(2:end)];
end

end
