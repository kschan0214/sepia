%% sepia_configuration_SWISMWI
%
% Description: a script contains all the SWI/SMWI method names and related function names
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 2 August 2022
% Date modified: 
%
% DO NOT change the variable name
% DO NOT change the order of the entities, add a new one at the end instead
% The order of the methods MUST be consistent across all variables
%
%% QSM dipole inversion methods available in SEPIA (MUST)
% Specify the names of the methods
methodSWISMWIName    	= {'SWI (2D Hamming)',...
                           'SMWI',...
                           };
                       % in future, add panel of new method here

%% Wrapper to access algorithms (MUST)
% Specify the name of the wrapper function to connect SEPIA and the
% algorithm
wrapper_SWISMWI_function	= {'Wrapper_SWISMWI_SWI_2DHamming',...
                               'Wrapper_SWISMWI_SMWI',...
                               } ;
                       % in future, add panel of new method here

%% function that create the method specific panel (Optional, for GUI)
% Specify the name of the function of the panel design 
function_SWISMWI_method_panel = {'sepia_handle_panel_swismwi_SWI_2DHamming',...
                                 'sepia_handle_panel_swismwi_SMWI',...
                                 };
                             % in future, add panel of new method here
                             
%% functions that generate and read algorithm parameters (Optional, for GUI)
% Specify the name of the function that export/import parameters to/from
% output config file
config_SWISMWI_function	= {'get_set_swismwi_swi_2dhamming',...
                           'get_set_swismwi_smwi',...
                           };
                       % in future, add panel of new method here
