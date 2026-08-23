function availableFileList = MaskWrapper(maskFullName, inputDir, sepia_header, algorParam, availableFileList, outputFileList, outputNiftiTemplate)
%MASKWRAPPER Summary of this function goes here
%   Detailed explanation goes here

sepia_universal_variables;

disp('---------------');
disp('Signal mask ');
disp('---------------');

isBET               = algorParam.general.isBET;
if isfield(algorParam.general, 'brain_extraction_method')
    brainExtractMethod  = algorParam.general.brain_extraction_method;
else
    brainExtractMethod = skullstrippingMethod{1};
end
if strcmp(brainExtractMethod,skullstrippingMethod{1})
    fractional_threshold    = algorParam.general.fractional_threshold;
    gradient_threshold      = algorParam.general.gradient_threshold;
end

matrixSize  = sepia_header.matrixSize;
voxelSize   = sepia_header.voxelSize;

mask        = [];
maskList    = dir(fullfile(inputDir,'*mask*nii*'));

% Scenario: No specified mask file + No check BET + there is a file called mask in the input directory
if isempty(maskFullName) && ~isempty(maskList) && ~isBET
    
    fprintf('No mask file is specified but a mask file is found in the input directory: %s\n',fullfile(inputDir, maskList(1).name));
    disp('Trying to load the file as signal mask');
    
    maskFullName = fullfile(inputDir, maskList(1).name);
end

% Scenario: User provided a mask file or above scenario was satified
if ~isempty(maskFullName)
    
    % load mask file
    mask = load_nii_img_only(maskFullName) > 0;
    
    % make sure the mask has the same dimension as other input data
    if ~isequal(size(mask),matrixSize)
        disp('The file does not have the same dimension as other images.')
        mask = [];
    else
        availableFileList.mask = maskFullName;
        disp('Mask file is checked.');
    end
end

% if no mask is found then display the following message
if isempty(mask) && ~isBET
    disp('No mask data is loaded. Using FSL BET to obtain brain mask.');
end
    
% if BET is checked or no mask is found, run FSL's bet
if isempty(mask) || isBET
    
    magn = load_nii_img_only(availableFileList.magnitude);
    mag_e1 = magn(:,:,:,1);

    % for synthstrip
    [temp_dir,~,~] = fileparts(outputFileList.maskBrain);
    temp_nii = fullfile(temp_dir,'temp.nii.gz');

    switch brainExtractMethod
        case skullstrippingMethod{1}    % MEDI implementation of BET
    
            sepia_addpath('MEDI');
            
            disp('Performing FSL BET...');
            % Here uses MEDI toolboxes MEX implementation
            mask = BET(mag_e1,matrixSize,voxelSize,fractional_threshold,gradient_threshold);
            disp('Signal mask is obtained.');

            fprintf('Saving signal mask...')
            save_nii_quick(outputNiftiTemplate,mask, outputFileList.maskBrain);

            fprintf('Done!\n');
            

        case skullstrippingMethod{3}    % synthstrip

            save_nii_quick(outputNiftiTemplate,mag_e1, temp_nii);

            cmd = sprintf('mri_synthstrip -i %s -m %s',temp_nii,outputFileList.maskBrain);

            status = system(cmd);
            if status ~= 0
                error('Failed running SynthStrip in the system. Please check if the tool is available in the PATH environment ot use other methods instead.');
            end
            delete(temp_nii);

        case skullstrippingMethod{4}    % synthstrip-no-CSF

            save_nii_quick(outputNiftiTemplate,mag_e1, temp_nii);

            cmd = sprintf('mri_synthstrip -i %s -m %s --no-csf',temp_nii,outputFileList.maskBrain);

            status = system(cmd);
            if status ~= 0
                error('Failed running SynthStrip in the system. Please check if the tool is available in the PATH environment ot use other methods instead.');
            end
            delete(temp_nii);
       
        case skullstrippingMethod{2}    % Otsu thresholding
            
            disp("Performing Otsu's Method Thresholding...");

            % Simple, single threshold segmentation
            % histogramCounts = histcounts(mag_e1);
            % level = otsu(histogramCounts);
            % mask = mag_e1 >= level;

            % Slightly more advanced (and robust) multilevel built-in
            thresh = multithresh(magn,3);
            labels = imquantize(mag_e1,thresh);

            % Select largest component of the labeled regions
            brain=labels==3;
            brain_props = regionprops(brain, 'PixelIdxList', 'Area');
            [~,bi] = max([brain_props.Area]);

            csf=labels==2;
            csf_props = regionprops(csf, 'PixelIdxList', 'Area');
            [~,ci] = max([csf_props.Area]);

            mask = zeros(size(mag_e1));
            mask(brain_props(bi).PixelIdxList) = 1;
            mask(csf_props(ci).PixelIdxList) = 1;



            fprintf('Saving signal mask...')
            save_nii_quick(outputNiftiTemplate,mask, outputFileList.maskBrain);

            fprintf('Done!\n');
    end

    if exist(outputFileList.maskBrain,'file')
        fprintf('Brain extraction Done!\n');
        availableFileList.mask = outputFileList.maskBrain;
    else
        error('No signal mask is found. QSM cannot be run without a signal mask.');
    end
end

end