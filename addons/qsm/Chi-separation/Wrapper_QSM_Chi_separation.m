%% [chi,chi_para,chi_dia] = Wrapper_QSM_Chi_separation(localField,mask,matrixSize,voxelSize,algorParam,headerAndExtraData)
%
% Input
% --------------
% localField    : local field map (tissue fields), in Hz
% mask          : signal mask
% matrixSize    : size of the input image
% voxelSize     : spatial resolution of each dimension of the data, in mm
% algorParam    : structure contains fields with algorithm-specific parameter(s)
% headerAndExtraData : structure contains extra header info/data for the algorithm
%
% Output
% --------------
% chi           : magnetic susceptibility map, in ppm
% chi_para      : Paramagnetic magnetic susceptibility map, in ppm
% chi_dia       : Diamagnetic magnetic susceptibility map, in ppm
%
% Description: This is a wrapper function to access Chi-Separation toolbox for SEPIA
% Required: importONNXNetwork requires the Deep Learning Toolbox Converter for ONNX Model Format support package. To install this support package, use the Add-On Explorer.)
%
% Date created: 02 Dec 2025 by Taechang Kim @ SNU (sakkar2@snu.ac.kr)
% Date modified: 22 August 2026 (KC)
%
% You can change the name of the function but DO NOT change the input/output variables
function [chi,chi_para,chi_dia] = Wrapper_QSM_Chi_separation(localField,mask,matrixSize,voxelSize,algorParam,headerAndExtraData)
% load some constants 
sepia_universal_variables;

% setup path
setup_Chi_sepnet_environment;
addpath(genpath(fullfile(home_directory,'functions')));
addpath(genpath(fullfile(home_directory,'models')));
addpath(genpath(fullfile(home_directory,'utils')));

% get algorithm parameters, if user doesn't specify them then set some default values
% algorParam = check_and_set_algorithm_default(algorParam);
solver         = algorParam.qsm.solver;
Dr             = algorParam.qsm.Dr;
if isfield(algorParam.qsm,'R2s')
    R2s_path = algorParam.qsm.R2s;
else
    R2s_path = ''; 
end
if isfield(algorParam.qsm,'R2')
    R2_path = algorParam.qsm.R2;
else
    R2_path = ''; 
end

% get extra data such as magnitude/weights/B0 direction/TE/etc.
headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);
params.b0_dir      = headerAndExtraData.sepia_header.B0_dir;
params.delta_TE    = headerAndExtraData.sepia_header.delta_TE;
params.TE          = headerAndExtraData.sepia_header.TE;
params.CF          = headerAndExtraData.sepia_header.CF;
params.B0_strength = headerAndExtraData.sepia_header.B0;
params.voxel_size  = voxelSize;
params.matrix_size = matrixSize;
params.Dr          = Dr;

iMag          = get_variable_from_headerAndExtraData(headerAndExtraData, 'magnitude', matrixSize);
iPhase        = get_variable_from_headerAndExtraData(headerAndExtraData, 'phase', matrixSize);
weights       = get_variable_from_headerAndExtraData(headerAndExtraData, 'weights', matrixSize);

%% Preparation

% R2prime
if strcmp(R2s_path,'')
    % R2* mapping: reuse the cached R2* map if one is already available,
    % otherwise compute it using ARLO
    if isfield(headerAndExtraData,'availableFileList') && isfield(headerAndExtraData.availableFileList,'r2s') && exist(headerAndExtraData.availableFileList.r2s,'file')
        disp('R2star path was not entered. R2* map is already available. Loading it from disk...');
        R2s = get_variable_from_headerAndExtraData(headerAndExtraData, 'r2s', matrixSize);
    else
        disp('R2star path was not entered. R2star is created using ARLO (MEDI toolbox)')
        sepia_addpath('MEDI');
        R2s = arlo(params.TE,iMag);
    end
else
    headerAndExtraData.availableFileList.r2s = R2s_path;
    R2s = get_variable_from_headerAndExtraData(headerAndExtraData, 'r2s', matrixSize);
end

if strcmp(R2_path,'')
    disp('R2 path was not entered. Chi-sepnet-R2* is performed, which only utilizes GRE data.')
    solver = 'Chi-sepnet-R2*';
    R2n = R2s .* mask;
    HaveR2prime = 0;
else
    % use get_variable_from_headerAndExtraData to get data that have the same matrix size as the loaded data 
    headerAndExtraData.availableFileList.R2 = R2_path;
    R2 = get_variable_from_headerAndExtraData(headerAndExtraData, 'R2', matrixSize);
    R2p = R2s - R2; R2p(R2p < 0) = 0;
    R2n = R2p .* mask;
    HaveR2prime = 1;
end

% CSF mask
sepia_addpath('MEDI'); % extract_CSF require MEDI toolbox
Mask_CSF = extract_CSF(R2s,mask,voxelSize)>0;

% weight function 
% 20260822 KC: in case a weight map is provided but not iPhase (e.g., using the QSM panel)
if ~isempty(iPhase) && ~isempty(iMag) % if phase exists then the priority is to use phase data as Taaecheng suggested
    [~, N_std] = Preprocessing4Phase(iMag, iPhase);
else
    N_std                   = 1./weights;   % invert weight to emulate N_std
    N_std(~isfinite(N_std)) = 0;
end

%% Display algorithm parameters + main
switch solver
   case 'Chi-separation-MEDI'
        disp('Chi-separation-MEDI is running');
        
        mag = sqrt(sum(iMag.^2,4)) .* mask;

        params.lambda = 1;
        params.lambda_CSF = 1;
        
        option_data.qsm = [];
        option_data.mask_CSF = Mask_CSF;
        option_data.N_std = N_std;
        option_data.wG = [];
        option_data.wG_r2p = [];
        option_data.mask_FastRelax = zeros(matrixSize);
        option_data.mask_SlowRelax = zeros(matrixSize);

        [x_para, x_dia, x_tot] = chi_sep_MEDI(mag, localField, R2n, N_std, mask, params, option_data);

    case 'Chi-separation-iLSQR'
        disp('Chi-separation-iLSQR is running');
        
        mag = sqrt(sum(iMag.^2,4)) .* mask;

        option_data.qsm = [];
        option_data.N_std = N_std;

        [x_para, x_dia, x_tot] = chi_sep_iLSQR(mag, localField, R2p, mask, params, option_data);

    case 'Chi-sepnet-R2*'
        disp('Chi-sepnet-R2* is running');

        % setup_Chi_sepnet_environment
        [x_para, x_dia, x_tot, ~, ~] = chi_sepnet_general_new_wResolGen(home_directory, localField, R2n, mask, params.Dr, ...
            params.b0_dir, params.CF, params.voxel_size, HaveR2prime, params.B0_strength, 0, 0.19, 'sinc', 15, 'hann');

    case 'Chi-sepnet-R2'''
        disp('Chi-sepnet-R2'' is running');

        % setup_Chi_sepnet_environment
        [x_para, x_dia, x_tot, ~, ~] = chi_sepnet_general_new_wResolGen(home_directory, localField, R2n, mask, params.Dr, ...
            params.b0_dir, params.CF, params.voxel_size, HaveR2prime, params.B0_strength, 0, 0.19, 'sinc', 15, 'hann');
        
 end

chi         = x_tot;
chi_para    = x_para;
chi_dia     = x_dia;

% 20260822 KC: move to I/O layer
% outputNiftiTemplate = load_untouch_nii(headerAndExtraData.availableFileList.localField);
% save_nii_quick(outputNiftiTemplate, x_para, fullfile(headerAndExtraData.outputDirectory,'Sepia_ChiPara.nii.gz'));
% save_nii_quick(outputNiftiTemplate, x_dia, fullfile(headerAndExtraData.outputDirectory,'Sepia_ChiDia.nii.gz'));

end
