%% addon_config
%
% Description: Addon script to incorporate CLEAR-SWI to the SEPIA framework
%
% Korbinian Eckstein @ UQ
% korbinian90@gmail.com
% Date created: 26 September 2022
% Date modified: 
%
% DO NOT change the variable names, update the variable value instead
%
%% main

% This name will be used thorough the SEPIA framework for usage
addons.method = 'CLEAR-SWI';

% Specify the filename of the wrapper function (without extension)
addons.wrapper_function	= 'Wrapper_CLEARSWI';

% Specify the filename of the GUI method panel function (without
% extension) 
% (optional): if no GUI support is provided, then set it to empty, i.e.
% addons.gui_method_panel	= [];
addons.gui_method_panel	= 'sepia_handle_panel_clearswi';

% Specify the filename of the function to export/import algorithm
% parameters (without extension) (optional)
% (optional): if no tunable algorithm parameter, then set it to empty, i.e.
% addons.config_function	= [];
addons.config_function	= 'get_set_clearswi';