%% addon_config
%
% Description: Addon script to incorporate new method to the SEPIA framework
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 27 June 2020
% Date modified: 
%
% DO NOT change the variable names, update the variable value instead
%
%% main

% This name will be used thorough the SEPIA framework for usage
addons.method = 'ROMEO';

% Specify the filename of the wrapper function (without extension)
addons.wrapper_function	= 'Wrapper_Unwrap_ROMEO_singleEcho';

% Specify if the exclusion of unreliable voxels option in the GUI is
% enabled. The method must return the ture phase value if it is enabled.
% Enable: 'on'
% Disable: 'off'
addons.gui_exclude_voxel= 'on';

% Specify the filename of the GUI method panel function (without
% extension) 
% (optional): if no GUI support is provided, then set it to empty, i.e.
% addons.gui_method_panel	= [];
addons.gui_method_panel	= 'sepia_handle_panel_Unwrap_ROMEO_singleEcho';

% Specify the filename of the function to export/import algorithm
% parameters (without extension) (optional)
% (optional): if no tunable algorithm parameter, then set it to empty, i.e.
% addons.config_function	= [];
addons.config_function	= 'get_set_Unwrap_ROMEO_singleEcho';
