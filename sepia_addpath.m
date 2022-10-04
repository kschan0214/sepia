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
% Date modified: 22 Feb 2022 (v1.0)
% Date modified: 3 August 2022 (v1.1)
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
    if ~exist('MEDI_HOME','var')
        MEDI_HOME = [];
    end
    if ~exist('STISuite_HOME','var')
        STISuite_HOME = [];
    end
    if ~exist('FANSI_HOME','var')
        FANSI_HOME = [];
    end
    if ~exist('SEGUE_HOME','var')
        SEGUE_HOME = [];
    end
    if ~exist('MRITOOLS_HOME','var')
        MRITOOLS_HOME = [];
    end
    if ~exist('MRISC_HOME','var')
        MRISC_HOME = [];
    end
        
    CheckPathValidity(MEDI_HOME,STISuite_HOME,FANSI_HOME,SEGUE_HOME,MRITOOLS_HOME,MRISC_HOME);
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
addpath(genpath(fullfile(SEPIA_HOME, 'addons')));

% misc directory accommodates some non-toolbox specific algorithms
misc_dir                = fullfile(SEPIA_HOME,  'misc');
misc_swi_smwi_dir       = fullfile(misc_dir,    'swi_smwi');
misc_r2s_dir            = fullfile(misc_dir,    'r2s_mapping');

addpath(genpath( misc_swi_smwi_dir));
addpath(genpath( misc_r2s_dir));

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
        
    case 'mritools'
        addpath(genpath(MRITOOLS_HOME));
        
    case 'mrisc'
        addpath(genpath(MRISC_HOME));

end

end

%% check the following paths exist or not
function CheckPathValidity(MEDI_HOME,STISuite_HOME,FANSI_HOME,SEGUE_HOME,MRITOOLS_HOME,MRISC_HOME)

if exist(MEDI_HOME,'dir')~=7
%     warning('Please specify a correct path for MEDI toolbox in SpecifyToolboxesDirectory.m');
    disp('The directory to MEDI toolbox does not exist.')
    warning('All functions related to MEDI toolbox cannot be used.');
end

if exist(STISuite_HOME,'dir')~=7
%     warning('Please specify a correct path for STI Suite in SpecifyToolboxesDirectory.m');
    disp('The directory to STI Suite does not exist.')
    warning('All functions related to STI Suite cannot be used.');
end

if exist(FANSI_HOME,'dir')~=7
%     warning('Please specify a correct path for FANSI toolbox in SpecifyToolboxesDirectory.m');
    disp('The directory to FANSI toolbox does not exist.')
    warning('All functions related to FANSI toolbox cannot be used.');
end

if exist(SEGUE_HOME,'dir')~=7
%     warning('Please specify a correct path for SEGUE in SpecifyToolboxesDirectory.m');
    disp('The directory to SEGUE toolbox does not exist.')
    warning('All functions related to SEGUE cannot be used.');
end
   
if exist(MRITOOLS_HOME,'dir')~=7
%     warning('Please specify a correct path for ROMEO in SpecifyToolboxesDirectory.m');
    disp('The directory to mritools does not exist.')
    warning('All functions related to ROMEO/CLEARSWI cannot be used.');
end

if exist(MRISC_HOME,'dir')~=7
%     warning('Please specify a correct path for MRI Susceptibility Calculation in SpecifyToolboxesDirectory.m');
    disp('The directory to MRI Susceptibility Calculation Methods does not exist.')
    warning('All functions related to MRI Susceptibility Calculation cannot be used.');
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