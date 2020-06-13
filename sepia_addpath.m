%% sepia_addpath(method)
%
% Description: remove all SEPIA related directories from PATH and add
% only the directory(ies) related to the input 'method'
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 10 April 2018
% Date modified: 24 August 2018
% Date modified: 5 June 2019
% Date modified: 9 June 2020 (v0.8.0)
%
function sepia_addpath(method, isStartCheck)

if nargin < 2
    isStartCheck = 0;
end
if nargin < 1
    method = 'None';
end

% specify the toolbox(es) directory
SpecifyToolboxesDirectory;
if isStartCheck
    CheckPathValidity(MEDI_HOME,STISuite_HOME,FANSI_HOME,SEGUE_HOME);
end

% get SEPIA_HOME from this file 
SEPIA_HOME = fileparts(mfilename('fullpath'));

% disable warning related to remove path
warning('off');
% remove all related paths except it root directory
rmpath(genpath(SEPIA_HOME));
addpath(SEPIA_HOME);
% enable warning
warning('on')

% essential path
addpath(fullfile(SEPIA_HOME,'configuration'));
addpath(genpath(fullfile(SEPIA_HOME, 'gui_func')));
addpath(genpath(fullfile(SEPIA_HOME, 'wrapper')));
addpath(genpath(fullfile(SEPIA_HOME, 'utils')));
addpath(genpath(fullfile(SEPIA_HOME, 'addon')));

% misc directory accommodates some non-toolbox specific algorithms
misc_dir                = fullfile(SEPIA_HOME,  'misc');
misc_swi_smwi_dir       = fullfile(misc_dir,    'swi_smwi');

addpath(misc_swi_smwi_dir);

% if method is given then add them to PATH
sepia_universal_variables;

switch lower(method)
        
    case 'medi'
        addpath(genpath(MEDI_HOME));
        
    case 'stisuite'
        addpath_STIsuitev3(STISuite_HOME);
        
    case 'segue'
        addpath(genpath(SEGUE_HOME));
        
    case 'fansi'
        addpath(genpath(FANSI_HOME));

end

end

%% check the following paths exist or not
function CheckPathValidity(MEDI_HOME,STISuite_HOME,FANSI_HOME,SEGUE_HOME)

if exist(MEDI_HOME,'dir')~=7
    warning('Please specify a correct path for MEDI toolbox in SpecifyToolboxesDirectory.m');
    warning('All functions related to MEDI toolbox cannot be used.');
end

if exist(STISuite_HOME,'dir')~=7
    warning('Please specify a correct path for STI Suite in SpecifyToolboxesDirectory.m');
    warning('All functions related to STI Suite cannot be used.');
end

if exist(FANSI_HOME,'dir')~=7
    warning('Please specify a correct path for FANSI toolbox in SpecifyToolboxesDirectory.m');
    warning('All functions related to FANSI toolbox cannot be used.');
end

if exist(SEGUE_HOME,'dir')~=7
    warning('Please specify a correct path for SEGUE in SpecifyToolboxesDirectory.m');
    warning('All functions related to SEGUE cannot be used.');
end
    
end

%% Special for STI suite
function addpath_STIsuitev3(STISuite_HOME)

addpath(fullfile(STISuite_HOME,'Core_Functions_P'));
addpath(fullfile(STISuite_HOME,'GUI_Functions_P'));
addpath(fullfile(STISuite_HOME,'Support_Functions'));
addpath(genpath(fullfile(STISuite_HOME,'Support_Functions','qsm_kiwi_1')));
addpath(genpath(fullfile(STISuite_HOME,'Support_Functions','SpaRSA')));
addpath(genpath(fullfile(STISuite_HOME,'Support_Functions','wavelet_src')));

end