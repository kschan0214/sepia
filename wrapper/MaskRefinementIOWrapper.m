%% availableFileList = MaskRefinementWrapper(Mask,algorParam,headerAndExtraData)
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

[mask_refined,r2s,residual] = MaskRefinementMacro(mask,algorParam,headerAndExtraData);

if not(isempty(r2s)) && isempty(availableFileList.r2s)
    disp('Saving R2* map.');
    save_nii_quick(outputNiftiTemplate, r2s, outputFileList.r2s);
    availableFileList.r2s = outputFileList.r2s;
end

if not(isempty(residual)) && isempty(availableFileList.relativeResidual)
    disp('Saving relative residual map.');
    save_nii_quick(outputNiftiTemplate, residual, outputFileList.relativeResidual);
    availableFileList.relativeResidual = outputFileList.relativeResidual;
end

disp('Saving refined brain mask.');
save_nii_quick(outputNiftiTemplate, mask_refined, outputFileList.maskReliable);
availableFileList.maskReliable = outputFileList.maskReliable;

end