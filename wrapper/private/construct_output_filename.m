%% outputFileList = construct_output_filename(outputDir, ouputPrefix)
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
function [outputFileList,ouputPrefix] = construct_output_filename(outputDir, ouputPrefix, algorParam)

% phase related
outputFileList.phaseRadian      = fullfile(outputDir, [ouputPrefix 'part-phase_rad.nii.gz']);
outputFileList.phaseReversed    = fullfile(outputDir, [ouputPrefix 'part-phase_reverse.nii.gz']);
outputFileList.phaseEddyCorr    = fullfile(outputDir, [ouputPrefix 'part-phase_bipolarcorr.nii.gz']);
outputFileList.unwrappedPhase   = fullfile(outputDir, [ouputPrefix 'part-phase_unwrapped.nii.gz']);

% standard output
outputFileList.totalField       = fullfile(outputDir, [ouputPrefix 'fieldmap.nii.gz']);
outputFileList.localField       = fullfile(outputDir, [ouputPrefix 'localfield.nii.gz']);
outputFileList.QSM              = fullfile(outputDir, [ouputPrefix 'Chimap.nii.gz']);

% use for regularisation
outputFileList.weights          = fullfile(outputDir, [ouputPrefix 'weights.nii.gz']);
outputFileList.fieldmapSD       = fullfile(outputDir, [ouputPrefix 'noisesd.nii.gz']);
outputFileList.relativeResidual	= fullfile(outputDir, [ouputPrefix 'relativeresidual.nii.gz']);
outputFileList.relativeResidualWeights	= fullfile(outputDir, [ouputPrefix 'relativeresidualweights.nii.gz']);

% derived masks
outputFileList.maskBrain        = fullfile(outputDir, [ouputPrefix 'mask_brain.nii.gz']);
outputFileList.maskReliable     = fullfile(outputDir, [ouputPrefix 'mask_reliable.nii.gz']);
outputFileList.maskLocalField 	= fullfile(outputDir, [ouputPrefix 'mask_localfield.nii.gz']);
outputFileList.maskQSM          = fullfile(outputDir, [ouputPrefix 'mask_QSM.nii.gz']);
outputFileList.maskRef      	= fullfile(outputDir, [ouputPrefix 'mask_referenceregion.nii.gz']);
outputFileList.maskRefine       = fullfile(outputDir, [ouputPrefix 'mask_refine.nii.gz']);

% R2*
outputFileList.r2s              = fullfile(outputDir, [ouputPrefix 'R2starmap.nii.gz']);
outputFileList.t2s              = fullfile(outputDir, [ouputPrefix 'T2starmap.nii.gz']);
outputFileList.s0               = fullfile(outputDir, [ouputPrefix 'S0map.nii.gz']);

% misc
outputFileList.phase_bipolar    = fullfile(outputDir, [ouputPrefix 'bipolar_phase.nii.gz']);

if algorParam.general.isDenoise
ouputPrefix = strcat(ouputPrefix,'denoised_');
outputFileList.magDenoise       = fullfile(outputDir, [ouputPrefix 'part-mag.nii.gz']);
outputFileList.phaseDenoise     = fullfile(outputDir, [ouputPrefix 'part-phase.nii.gz']);
outputFileList.sigma            = fullfile(outputDir, [ouputPrefix 'sigma.nii.gz']);
outputFileList.snrgain          = fullfile(outputDir, [ouputPrefix 'SNRgain.nii.gz']);
outputFileList.P                = fullfile(outputDir, [ouputPrefix 'P.nii.gz']);
end

if algorParam.general.isUpsample
ouputPrefix = strcat(ouputPrefix,'upsampled_');
outputFileList.magUpsample      = fullfile(outputDir, [ouputPrefix 'part-mag.nii.gz']);
outputFileList.phaseUpsample    = fullfile(outputDir, [ouputPrefix 'part-phase.nii.gz']);
outputFileList.maskUpsample     = fullfile(outputDir, [ouputPrefix 'mask_upsampled.nii.gz']);
outputFileList.sepiaHeaderUpsample = fullfile(outputDir, [ouputPrefix 'sepia_header.mat']);
end

end