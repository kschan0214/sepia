%% SpecifyAtlasDirectory.m
%
% Description: Specify directories of the atlas
%
% Kwok-Shing Chan @ MGH
% kchan2@mgh.harvard.edu
% Date created: 9 October 2023
% Date modified: 
%
%% Specify the directories of the atlas here
% please do not change
SEPIA_ANALYSIS_SEGMENTATION_dir = fullfile(SEPIA_HOME,'analysis','segmentation');
SEPIA_ATLAS_dir                 = fullfile(SEPIA_HOME,'atlas');

% These are the default directories if you use the 'download_atlas.sh' shell script
% If you prefer the atlases to be stored at a different paths, you may
% modified the paths here

% AHEAD atlas
AHEAD_ATLAS_HOME                = fullfile(SEPIA_ATLAS_dir,'AHEAD_atlas');
% CIT168 atlas
CIT168_reinf_learn_ATLAS_HOME   = fullfile(SEPIA_ATLAS_dir,'CIT168_Reinf_Learn_v1.1.0');
% MUSUS-100 atlas
MuSus100_ATLAS_HOME             = fullfile(SEPIA_ATLAS_dir,'MuSus-100_Atlas');