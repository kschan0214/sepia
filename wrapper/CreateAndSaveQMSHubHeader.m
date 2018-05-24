%% CreateAndSaveQMSHubHeader(input,outputDir,b0_input,b0dir_input,voxelSize_input,te_input,flag)
%
% Input
% --------------
% input             : either DICOM directory or NIfTI file fullpath
% outputDir         : directory that store the output header file
% b0_input          : user field strength input
% b0dir_input       : user field direction input
% voxelSize_input   : user voxel size input
% te_input          : user echo times input
% flag              : either 'nifti' or 'dicom', indicator of input type
%
% Description: Utility function to save header required by qsm_hub from
% DICOM, NIfTi or user input
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 24 May 2018
% Date last modified:
%
%
function CreateAndSaveQMSHubHeader(input,outputDir,b0_input,b0dir_input,voxelSize_input,te_input,flag)
gyro = 42.57747892;

switch flag
    case 'dicom'
    % for DICOM input, no user input will be saved
        [~,voxelSize,matrixSize,CF,delta_TE,TE,B0_dir]=Read_DICOM(input);
        
        B0 = CF/(gyro*1e6);
        
    case 'nifti'
    % for NIfTI input, header will be read from the file
    % user input is allowed
        inputNifti = load_untouch_nii(input);
        
        % create synthetic header in case no qsm_hub's header is found
        [B0,B0_dir,voxelSize,matrixSize,TE,delta_TE,CF]=SyntheticQSMHubHeader(inputNifti);
        
        % check user input
        if ~isempty(b0_input)
            B0 = b0_input;
            CF = B0 * gyro * 1e6;
        end
        if ~isempty(b0dir_input)
            b0dir_input = b0dir_input./norm(b0dir_input);
            B0_dir = b0dir_input;
        end
        if ~isempty(voxelSize_input)
            voxelSize = voxelSize_input;
        end
        if ~isempty(TE)
            TE = te_input;
            delta_TE = TE(2)-TE(1);
        end
        
end

save([outputDir filesep 'squirrel_header.mat'],'voxelSize','matrixSize','CF','delta_TE',...
        'TE','B0_dir','B0');

end

