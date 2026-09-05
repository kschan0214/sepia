%% availableFileList = MaskRefinementWrapper(Mask,algorParam,headerAndExtraData)
%
% Input
% --------------
% Mask              : Logical brain mask to be refined
% RefinementMap     : Map or data to use for refinement. This will depend
%                     on the refinement strategy, and can be of the type:
%                       - R2star map
%                       - Local fieldmap
%                       - Noise map
% algorParam        : Struct containing the method and method specific parameters
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
function availableFileList = MaskRefinementWrapper(sepia_header, algorParam, availableFileList, outputFileList, outputNiftiTemplate)

sepia_universal_variables;

TE           = sepia_header.TE;
voxelSize    = sepia_header.voxelSize;
refineMethod = algorParam.msk.refineMethod;
threshold    = algorParam.msk.threshold;

disp('---------------');
disp('Mask refinement');
disp('---------------');
fprintf('Refining brain mask...');

mask        = double(load_nii_img_only(availableFileList.mask));

switch lower(refineMethod)

    case {'r2s-refine'} % Original mask refinement code
        
        if isscalar(TE) 
            warning('Mask refinement using R2* only works for Multi-Echo data')
            return
        end

        disp('Refine brain using R2* info');
        % R2* map only needs to be computed once; reuse it if it is already
        % available, otherwise compute it and make it available for later use
        if isfield(availableFileList,'r2s') && exist(availableFileList.r2s,'file')
            disp('R2* map is already available. Loading it from disk...');
            r2s = double(load_nii_img_only(availableFileList.r2s));
        else
            magn = double(load_nii_img_only(availableFileList.magnitude));
            r2s  = R2star_trapezoidal(magn, TE);
    
            fprintf('Saving R2* map...');
            save_nii_quick(outputNiftiTemplate, r2s, outputFileList.r2s);
            fprintf('Done!\n');
    
            availableFileList.r2s = outputFileList.r2s;
        end
        
        maskRefined = refine_brain_mask_using_r2s(r2s,mask,voxelSize);
        
    case {lower(methodTwoPassName{2}), 'monoexponential','decay','mdm'}

        magn = load_nii_img_only(availableFileList.magnitude);
        fieldMap = load_nii_img_only(availableFileList.phase);
        totalField = load_nii_img_only(availableFileList.totalField);

        if availableFileList.r2s
            r2s = load_nii_img_only(availableFileList.r2s);
        else
            % multi-echo data
            r2s = R2star_trapezoidal(magn,TE); 
            % or arlo(headerAndExtraData.sepia_header.TE, magn);?
            save_nii_quick(outputNiftiTemplate,r2s, outputFileList.r2s);
            availableFileList.r2s = outputFileList.r2s;
        end
        relativeResidual   = ComputeResidualGivenR2sFieldmap(TE,r2s,totalField,magn.*exp(1i*fieldMap));
        maskRefined        = relativeResidual < threshold;
        %
        exclude_threshold = 0.5;
        relativeResidualWeights = relativeResidual;
        % clipping
        relativeResidualWeights(relativeResidualWeights>exclude_threshold) = exclude_threshold;
        % weightsRelativeResidual should be between [0,1]
        relativeResidualWeights = (exclude_threshold - relativeResidualWeights) ./ exclude_threshold;

        save_nii_quick(outputNiftiTemplate,relativeResidual,       outputFileList.relativeResidual);
        save_nii_quick(outputNiftiTemplate,relativeResidualWeights,outputFileList.relativeResidualWeights);
        fprintf('Done.\n');

        clear relativeResidual

        availableFileList.relativeResidual          = outputFileList.relativeResidual;
        availableFileList.relativeResidualWeights   = outputFileList.relativeResidualWeights;

    case {lower(methodTwoPassName{3}), 'mgf','magnitude gradient field'}
        if availableFileList.localField
            localField = load_nii_img_only(availableFileList.localField);
        else
            error('Magnitude Gradient Field masking is only possible when a fieldmap is available.')
        end
        maskRefined = GradientBasedThreshold(localField, mask, threshold);

        save_nii_quick(outputNiftiTemplate,maskRefined, outputFileList.maskReliable);
        availableFileList.maskReliable                 = outputFileList.maskReliable;

    case {lower(methodTwoPassName{4}),'noisemap','nstd'}

        if not(isempty(availableFileList.fieldmapSD))
            fieldmapSD  = load_nii_img_only(availableFileList.fieldmapSD);

            maskRefined = erode3d(mask, fieldmapSD);
        else
            warning('Cannot perform noise based mask refinement when no noisemap is provided')
        end

    otherwise
        warning('Invalid mask refinement strategy given, keeping mask unchanged.');
        maskRefined = mask;
end

fprintf('Saving refined brain mask...');
save_nii_quick(outputNiftiTemplate,maskRefined, outputFileList.maskReliable);
fprintf('Done.\n');

% update availableFileList
availableFileList.maskReliable = outputFileList.maskReliable;

end