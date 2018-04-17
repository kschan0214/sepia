%% qsm_hub_AddMethodPath(method)
%
% Description: remove all qsm_hub related directories from PATH and add
% only the directory(ies) related to the input 'method'
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 10 April 2018
% Date last modified: 
%
function qsm_hub_AddMethodPath(method)
% specify the qsm algorithm version
% TODO: may allow user to choose software versions
MEDI_version = 'MEDI_20171106';
STISuite_version = 'STISuitev3'; % STISuite_version = 'STISuitev2_2';

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
addpath([currDir filesep 'gui_func/']);
addpath([currDir filesep 'macro/']);

% here are three basic directories of qsm_hub
phaseUnwrapDir = [currDir filesep 'phase_unwrap/'];
bkgRemovalDir = [currDir filesep 'background_removal/'];
qsmAlgorithmDir = [currDir filesep 'qsm_algorithm/'];
addpath(phaseUnwrapDir);
addpath(bkgRemovalDir);
addpath(qsmAlgorithmDir);
addpath([phaseUnwrapDir 'utils']);
addpath([bkgRemovalDir 'utils']);
addpath([qsmAlgorithmDir 'utils']);


% these directories contains codes that are not algorithm-specific
utilsDir = [currDir filesep 'utils/'];
addpath(utilsDir);
addpath([utilsDir 'nifti/NIfTI_20140122/']);
addpath([utilsDir 'nifti/utils/']);
addpath([utilsDir 'ReadDICOMMEDI/' MEDI_version]);

% if method is given then add them to PATH
if nargin > 0
    switch lower(method)
        case 'bet'
            addpath(genpath([utilsDir MEDI_version]));
        case 'laplacian'
            addpath([phaseUnwrapDir 'unwrapLaplacian' filesep MEDI_version]);
        case 'laplacian_stisuite'
            addpath([phaseUnwrapDir 'unwrapLaplacian' filesep STISuite_version]);
            addpath(genpath([utilsDir STISuite_version]));
        case 'rg'
            addpath([phaseUnwrapDir 'unwrapRegionGrowing' filesep MEDI_version]);
        case 'gc'    
            addpath(genpath([phaseUnwrapDir 'unwrapGraphcut' filesep MEDI_version]));
        case 'jena'
            addpath([phaseUnwrapDir 'unwrapJena']);
        case 'lbv'
            addpath([bkgRemovalDir 'LBV' filesep MEDI_version]);
        case 'pdf'
            addpath([bkgRemovalDir 'PDF' filesep MEDI_version]);
            addpath(genpath([utilsDir MEDI_version]));
        case 'sharp'
            addpath([bkgRemovalDir 'SHARP']);
        case 'resharp'
            addpath([bkgRemovalDir 'RESHARP']);
        case 'vsharpsti'
            addpath([bkgRemovalDir 'VSHARP' filesep STISuite_version]);
            addpath(genpath([utilsDir STISuite_version]));
        case 'vsharp'
            addpath([bkgRemovalDir 'VSHARP' filesep 'vsharp_chan']);
        case 'iharperella'
            addpath([bkgRemovalDir 'iHARPERELLA' filesep STISuite_version]);
            addpath(genpath([utilsDir STISuite_version]));
        case 'tkd'
            addpath([qsmAlgorithmDir 'TKD']);
        case 'closedforml2'
            addpath([qsmAlgorithmDir 'closedFormL2']);
        case 'ilsqr'
            addpath([qsmAlgorithmDir 'iLSQR' filesep 'iLSQR_chan']);
        case 'stisuiteilsqr'
            addpath([qsmAlgorithmDir 'iLSQR' filesep STISuite_version]);
            addpath(genpath([utilsDir STISuite_version]));
        case 'fansi'
            addpath([qsmAlgorithmDir 'FANSI']);
        case 'ssvsharp'
            addpath([qsmAlgorithmDir 'SingleStep']);
        case 'star'
            addpath([qsmAlgorithmDir 'Star']);
            addpath(genpath([utilsDir STISuite_version]));
        case 'medi_l1'
            addpath([qsmAlgorithmDir 'MEDI' filesep MEDI_version]);
            addpath(genpath([utilsDir MEDI_version]));
    end
end

end