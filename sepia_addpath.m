%% sepia_addpath(method)
%
% Description: remove all qsm_hub related directories from PATH and add
% only the directory(ies) related to the input 'method'
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 10 April 2018
% Date modified: 24 August 2018
% Date modified: 5 June 2019
%
function sepia_addpath(method)
sepia_universal_variables;
methodUnwrapName = lower(methodUnwrapName);
methodBFRName    = lower(methodBFRName);
methodQSMName    = lower(methodQSMName);

% specify the toolbox(es) directory
SpecifyToolboxesDirectory;
CheckPathValidity(MEDI_dir,STISuite_dir,FANSI_dir,SEGUE_dir);

% get the full path of this file
fullName = mfilename('fullpath');
currDir = fileparts(fullName);

% disable warning related to remove path
warning('off');
% remove all related paths except it root directory
rmpath(genpath(currDir));
addpath(currDir);
% enable warning
warning('on')

% keep GUI and macro directories
addpath(genpath([currDir filesep 'gui_func/']));
addpath(genpath([currDir filesep 'wrapper/']));

% misc directory accommodates some non-toolbox specific algorithms
misc_dir = [currDir filesep 'misc/'];
misc_qsm_dir = [misc_dir filesep 'qsm_algorithm/'];
misc_bkgRemoval_dir = [misc_dir filesep 'background_removal/'];
misc_phaseUnwrap_dir = [misc_dir filesep 'phase_unwrap/'];
misc_swi_smwi_dir = [misc_dir filesep 'swi_smwi/'];
addpath(misc_swi_smwi_dir);

% these directories contains codes that are not algorithm-specific
utilsDir = [currDir filesep 'utils/'];
addpath(utilsDir);
addpath([utilsDir 'nifti/NIfTI_20140122/']);
% addpath([utilsDir 'nifti/utils/']);
addpath([utilsDir 'nifti/']);
addpath([utilsDir 'nifti/quaternions/']);

% if method is given then add them to PATH
if nargin > 0
    switch lower(method)
        % I/O
        case 'dicom'
            addpath(genpath(MEDI_dir));
        case 'bet'
            addpath(genpath(MEDI_dir));

        % phase unwrap
        case 'nonlinearfit'
            addpath(genpath(MEDI_dir));
            
        case methodUnwrapName{1} % 'laplacian (medi)'
            addpath(genpath(MEDI_dir));

        case methodUnwrapName{2} % 'laplacian (sti suite)'
            addpath(genpath(MEDI_dir));
            add_path_STIsuitev3(STISuite_dir);

        case methodUnwrapName{4} % 'region growing (medi)'
            addpath(genpath(MEDI_dir));

        case methodUnwrapName{5} % 'graphcut'    
            addpath(genpath(MEDI_dir));

        case methodUnwrapName{3} % '3d best path'
            addpath(genpath(MEDI_dir));
            addpath([misc_phaseUnwrap_dir filesep 'unwrapBestpath3D/']);
            
        case methodUnwrapName{6} % 'segue'
            addpath(genpath(SEGUE_dir));

        % background field removal
        case methodBFRName{1} % 'lbv'
            addpath(genpath(MEDI_dir));

        case methodBFRName{2} % 'pdf'
            addpath(genpath(MEDI_dir));

        case methodBFRName{4} % 'sharp'
            addpath([misc_bkgRemoval_dir filesep 'SHARP']);
            addpath(genpath(MEDI_dir));

        case methodBFRName{3} % 'resharp'
            addpath([misc_bkgRemoval_dir filesep 'RESHARP']);
            addpath(genpath(MEDI_dir));

        case methodBFRName{5} % 'vsharpstisuite'
            add_path_STIsuitev3(STISuite_dir);

        case methodBFRName{6} % 'vsharp'
            addpath([misc_bkgRemoval_dir filesep 'VSHARP_sepia']);

        case methodBFRName{7} % 'iharperella'
            add_path_STIsuitev3(STISuite_dir);

        % QSM
        case methodQSMName{1} % 'tkd'
            addpath([misc_qsm_dir filesep 'TKD']);

        case methodQSMName{2} % 'cfl2'
            addpath([misc_qsm_dir filesep 'closedFormL2']);
            
        case methodQSMName{3} % 'ndi'
            addpath([misc_qsm_dir filesep 'NDI']);

        case methodQSMName{5} % 'ilsqr'
            addpath([misc_qsm_dir filesep 'closedFormL2']);
            addpath([misc_qsm_dir filesep 'iLSQR_qsmhub']);

        case methodQSMName{4} % 'stisuiteilsqr'
            add_path_STIsuitev3(STISuite_dir);

        case methodQSMName{6} % 'fansi'
            addpath([misc_qsm_dir filesep 'FANSI']);
            addpath(genpath(FANSI_dir));

        case methodQSMName{7} % 'star'
            add_path_STIsuitev3(STISuite_dir);

        case methodQSMName{8} % 'medi_l1'
            addpath([misc_qsm_dir filesep 'MEDI_L1']);
            addpath(genpath(MEDI_dir));

    end
end

end

%% check the following paths exist or not
function CheckPathValidity(MEDI_dir,STISuite_dir,FANSI_dir,SEGUE_dir)

if exist(MEDI_dir,'dir')~=7
    warning('Please specify a correct path for MEDI toolbox in SpecifyToolboxesDirectory.m');
    warning('All functions related to MEDI toolbox cannot be used.');
end

if exist(STISuite_dir,'dir')~=7
    warning('Please specify a correct path for STI Suite in SpecifyToolboxesDirectory.m');
    warning('All functions related to STI Suite cannot be used.');
end

if exist(FANSI_dir,'dir')~=7
    warning('Please specify a correct path for FANSI toolbox in SpecifyToolboxesDirectory.m');
    warning('All functions related to FANSI toolbox cannot be used.');
end

if exist(SEGUE_dir,'dir')~=7
    warning('Please specify a correct path for SEGUE in SpecifyToolboxesDirectory.m');
    warning('All functions related to SEGUE cannot be used.');
end
    
end