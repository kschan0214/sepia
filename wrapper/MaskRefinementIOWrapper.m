%% availableFileList = MaskRefinementIOWrapper(sepia_header, algorParam, availableFileList, outputFileList, outputNiftiTemplate)
%
% Input
% --------------
% sepia_header              
% algorParam            : Struct containing the method and method specific parameters
% availableFileList                  
% outputFileList
% outputNiftiTemplate
%
% algorParam        
%
% Output
% --------------
% availableFileList      : Updated file-list
%
% Description: exclude unreliable mask voxels, based on various refinement 
% strategies. Possible strategies are:
%
%       - Monoexponential decay model   (requires R2star map)
%       - Magnitude Gradient Field      (requires local fieldmap)
%       - Noise map                     (requires noise map)
%
% Created by Patrick Fuchs
% Based on code by Kwok-shing Chan, Anita Karsa and Oliver C. Kiersnowski
% patrick.fuchs@uantwerpen.be
% Date created: 5 August 2025
% Date modified: 3 September 2026
%
function availableFileList = MaskRefinementIOWrapper(sepia_header, algorParam, availableFileList, outputFileList, outputNiftiTemplate)

sepia_universal_variables;

headerAndExtraData.sepia_header = sepia_header;
headerAndExtraData.availableFileList = availableFileList;
mask = double(load_nii_img_only(availableFileList.mask));

% derive the output directory from one of the (always present) output filenames
[outputDir,~,~] = fileparts(outputFileList.maskReliable);

[mask_refined,r2s,residual,gradientMagnitude,gradientStats] = MaskRefinementMacro(mask,algorParam,headerAndExtraData);

if not(isempty(r2s)) && (~isfield(availableFileList,'r2s') || isempty(availableFileList.r2s) )
    disp('Saving R2* map.');
    save_nii_quick(outputNiftiTemplate, r2s, outputFileList.r2s);
    save_json_sidecar(outputFileList.r2s, struct( ...
        'Description', 'R2* map estimated from multi-echo magnitude data for mask refinement.', ...
        'Units',       '1/s', ...
        'Sources',     {{get_relative_source_path(outputDir, availableFileList.magnitude)}}, ...
        'Method',       'Trapezoidal approximation'));
    availableFileList.r2s = outputFileList.r2s;
end

if not(isempty(residual)) && (~isfield(availableFileList,'relativeResidual') || isempty(availableFileList.relativeResidual) )
    disp('Saving relative residual map.');
    save_nii_quick(outputNiftiTemplate, residual, outputFileList.relativeResidual);
    save_json_sidecar(outputFileList.relativeResidual, struct( ...
        'Description', 'Relative residual between the measured and mono-exponential modelled magnitude decay, used for mask refinement.', ...
        'Units',       'ratio', ...
        'Sources',     {{get_relative_source_path(outputDir, availableFileList.magnitude)}}, ...
        'Parameters',  algorParam.msk));
    availableFileList.relativeResidual = outputFileList.relativeResidual;
end

if not(isempty(gradientMagnitude)) && (~isfield(availableFileList,'gradientMagnitude') || isempty(availableFileList.gradientMagnitude) )
    disp('Saving gradient magnitude map.');
    save_nii_quick(outputNiftiTemplate, gradientMagnitude, outputFileList.gradientMagnitude);
    save_json_sidecar(outputFileList.gradientMagnitude, struct( ...
        'Description', 'Magnitude of the local field map gradient, used by the Magnitude Gradient Field mask refinement strategy.', ...
        'Units',       'Hz/voxel', ...
        'Sources',     {{get_relative_source_path(outputDir, availableFileList.localField)}}, ...
        'Parameters',  algorParam.msk));
    availableFileList.gradientMagnitude = outputFileList.gradientMagnitude;
end

disp('Saving refined brain mask.');
save_nii_quick(outputNiftiTemplate, mask_refined, outputFileList.maskReliable);
maskReliableInfo = struct( ...
        'Description', 'Refined signal mask excluding unreliable voxels.', ...
        'Units',       'binary', ...
        'Sources',     {{get_relative_source_path(outputDir, availableFileList.mask)}}, ...
        'Parameters',  algorParam.msk);
if not(isempty(gradientStats))
    % Magnitude Gradient Field: record the mean/std/threshold actually used,
    % so the threshold applied to this specific mask can be traced back
    % without re-deriving it from the gradient magnitude map.
    maskReliableInfo.GradientStatistics = gradientStats;
end
save_json_sidecar(outputFileList.maskReliable, maskReliableInfo);
availableFileList.maskReliable = outputFileList.maskReliable;

end