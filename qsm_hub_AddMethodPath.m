%% qsm_hub_AddMethodPath(method)
%
% Description: remove all qsm_hub related directories from PATH and add
% only the directory(ies) related to the input 'method'
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 10 April 2018
% Date last modified: 24 August 2018
%
function qsm_hub_AddMethodPath(method)
% specify the toolbox(es) directory
MEDI_dir = '/home/mrphys/kwocha/Tools/squirrel/qsm_hub/MEDI_toolbox/MEDI_toolbox_20180625/';
STISuite_dir = '/home/mrphys/kwocha/Tools/squirrel/qsm_hub/STI_Suite/STISuite_V3.0/';
FANSI_dir = '/home/mrphys/kwocha/Tools/squirrel/qsm_hub/FANSI_toolbox/FANSI-toolbox-d33759b970790cc8754adc9d0398cc3d07546074/';


% get the full path of this file
fullName = mfilename('fullpath');
currDir = fileparts(fullName);

% disable warning related to remove path
warning('off');
% remove all qsm_hub related paths except it root directory
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

% these directories contains codes that are not algorithm-specific
utilsDir = [currDir filesep 'utils/'];
addpath(utilsDir);
addpath([utilsDir 'nifti/NIfTI_20140122/']);
addpath([utilsDir 'nifti/utils/']);
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
        case 'laplacian'
            addpath(genpath(MEDI_dir));

        case 'laplacian_stisuite'
            addpath(genpath(MEDI_dir));
            addpath(genpath(STISuite_dir));

        case 'regiongrowing'
            addpath(genpath(MEDI_dir));

        case 'graphcut'    
            addpath(genpath(MEDI_dir));

        case 'bestpath3d'
            addpath(genpath(MEDI_dir));
            addpath([misc_phaseUnwrap_dir filesep 'unwrapBestpath3D/']);

        % background field removal
        case 'lbv'
            addpath(genpath(MEDI_dir));

        case 'pdf'
            addpath(genpath(MEDI_dir));

        case 'sharp'
            addpath([misc_bkgRemoval_dir filesep 'SHARP']);
            addpath(genpath(MEDI_dir));

        case 'resharp'
            addpath([misc_bkgRemoval_dir filesep 'RESHARP']);
            addpath(genpath(MEDI_dir));

        case 'vsharpsti'
            addpath(genpath(STISuite_dir));

        case 'vsharp'
            addpath([misc_bkgRemoval_dir filesep 'VSHARP_qsmhub']);

        case 'iharperella'
            addpath(genpath(STISuite_dir));

        % QSM
        case 'tkd'
            addpath([misc_qsm_dir filesep 'TKD']);

        case 'closedforml2'
            addpath([misc_qsm_dir filesep 'closedFormL2']);

        case 'ilsqr'
            addpath([misc_qsm_dir filesep 'closedFormL2']);
            addpath([misc_qsm_dir filesep 'iLSQR_qsmhub']);

        case 'stisuiteilsqr'
            addpath(genpath(STISuite_dir));

        case 'fansi'
            addpath(genpath(FANSI_dir));

        case 'star'
            addpath(genpath(STISuite_dir));

        case 'medi_l1'
            addpath([misc_qsm_dir filesep 'MEDI_L1']);
            addpath(genpath(MEDI_dir));

    end
end

end