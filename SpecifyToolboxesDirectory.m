%% SpecifyToolboxesDirectory.m
%
% Description: Specify the toolboxes's directory
% The specificed path must contain the following directories as original one:
% MEDI toolbox  - functions/
% STI Suite     - Core_Functions_P/
% 
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 24 August 2018
% Date last modified:
%
%
%% Specify the directories of the toolbox here
fullName = mfilename('fullpath');
currDir = fileparts(fullName);

MEDI_dir = [currDir filesep 'MEDI_toolbox/MEDI_toolbox_20180625/'];
STISuite_dir = [currDir filesep 'STI_Suite/STISuite_V3.0/'];
FANSI_dir = [currDir filesep 'FANSI_toolbox/FANSI-toolbox-d33759b970790cc8754adc9d0398cc3d07546074/'];
