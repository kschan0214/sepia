%% availableFileList = MaskRefinementWrapper(Mask,algorParam,headerAndExtraData)
%
% Input
% --------------
% Mask               : Logical brain mask to be refined
% algorParam         : Struct containing the method and method specific parameters
% headerAndExtraData : Struct containting the maps and data necessary for
%                      mask refinement
%
% Output
% --------------
% RefinedMask       : Refined mask
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
% Date modified: 
%
function [mask_refined,r2s,residual] = MaskRefinementMacro(mask,algorParam,headerAndExtraData)

sepia_universal_variables;

voxelSize    = headerAndExtraData.sepia_header.voxelSize;
TE           = headerAndExtraData.sepia_header.TE;
matrixSize   = size(mask);

refineMethod = algorParam.msk.refineMethod;
threshold    = algorParam.msk.threshold;

disp('--------------------');
disp('Mask refinement step');
disp('--------------------');

mask_refined = mask;
r2s = [];
residual = [];

switch lower(refineMethod)

    case {'r2s-refine'} % Original mask refinement code
        disp('Refine brain using R2* information.');
        % R2* map only needs to be computed once; reuse it if it is already
        % available, otherwise compute it and make it available for later use
        if isscalar(TE) 
            warning('Mask refinement using R2* only works for Multi-Echo data')
            return
        end

        if( isempty(headerAndExtraData.availableFileList.magnitude))
            warning('Please specify a magnitude data (at least 3 echoes) if you want to use r2s based mask refinement.');
            warning('No mask refinement is done in this instance.');
            return
        else
            sepia_addpath('MEDI');
            magn        = get_variable_from_headerAndExtraData(headerAndExtraData, 'magnitude');
            if size(magn,4) < 3
                warning('Please specify a magnitude data (at least 3 echoes) if you want to use CSF as reference.');
                warning('No normalisation will be done on the susceptibility map in this instance.');
                return                
            else
                if isfield(headerAndExtraData.availableFileList,'r2s') && ...
                    exist(headerAndExtraData.availableFileList.r2s,'file')
                    disp('R2* map is already available. Loading it from disk...');
                    r2s = get_variable_from_headerAndExtraData(headerAndExtraData, 'r2s');
                else
                    r2s  = R2star_trapezoidal(magn, TE);
                end
                clear magn
            end
        end

        mask_refined = refine_brain_mask_using_r2s(r2s,mask,voxelSize);
        
    case {lower(methodTwoPassName{2}), 'monoexponential','decay','mdm'}
        disp('Refine brain using mono-exponential decay model.');
        magn = check_and_load(headerAndExtraData,'magnitude');
        phase = check_and_load(headerAndExtraData, 'phase');
        totalField = check_and_load(headerAndExtraData, 'totalField');

        if isfield(headerAndExtraData.availableFileList,'r2s') && ...
            exist(headerAndExtraData.availableFileList.r2s,'file')
            disp('R2* map is already available. Loading it from disk...');
            r2s = double(load_nii_img_only(headerAndExtraData.availableFileList.r2s));
        else
            r2s  = R2star_trapezoidal(magn, TE);
        end
        relativeResidual   = ComputeResidualGivenR2sFieldmap(TE,r2s,totalField,magn.*exp(1i*phase));
        mask_refined        = (relativeResidual < threshold) & mask;

    case {lower(methodTwoPassName{3}), 'mgf','magnitude gradient field'}
        disp('Refine brain using the magnitude of the gradient of the fieldmap.');
        fprintf(['Please cite:\nhttps://archive.ismrm.org/2024/3674.html for',...
                     ' MGF masking.\n'] )
        localField = check_and_load(headerAndExtraData, 'localField');
        if isempty(localField)
            warning('Magnitude Gradient Field masking is only possible when a fieldmap is available.')
            return
        end

        mask_refined = GradientBasedThreshold(localField, mask, threshold);

    case {lower(methodTwoPassName{4}),'noisemap','nstd'}
        disp('Refine brain using the noise map.');
        fieldmapSD = check_and_load(headerAndExtraData, 'fieldmapSD');
        if isempty(fieldmapSD)
            warning('Cannot perform noise based mask refinement when no noisemap is provided')
            return
        end
        mask_refined = erode3d(mask, fieldmapSD);

    otherwise
        warning('Invalid mask refinement strategy given, keeping mask unchanged.');
        return

end

end

function map = check_and_load(headerAndExtraData,variableName)
        % Check availability required inputs
        if isfield(headerAndExtraData.availableFileList,variableName) && ...
            exist(headerAndExtraData.availableFileList.(variableName),'file')
            map = double(load_nii_img_only(headerAndExtraData.availableFileList.(variableName)));
        else
            warning('Please specify a %s data if you want to use mask refinement.',variableName);
            warning('No mask refinement is done in this instance.');
            map = [];
            return
        end

end