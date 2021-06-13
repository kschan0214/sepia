%% sepia_configuration_EchoCombine
%
% Description: a script contains all the echo combine method names and related function names
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 5 June 2021 (v1.0)
% Date modified: 
%
% DO NOT change the variable name
% DO NOT change the order of the entities, add a new one at the end instead
% The order of the methods MUST be consistent across all variables
%
%% Background field removal method available in SEPIA (MUST)
% Specify the names of the methods
methodEchoCombineName	= {'Optimum weights',...
                           'MEDI nonlinear fit',...
                           };
%                            'MEDI nonlinear fit (Bipolar, testing)'};
                       % in future, add panel of new method here

%% Wrapper to access EchoCombine algorithms (MUST)
% Specify the name of the wrapper function to connect SEPIA and the
% algorithm
wrapper_EchoCombine_function	= {'Wrapper_EchoCombine_OptimumWeight',...
                                   'Wrapper_EchoCombine_MEDInonlinearfit',...
                                   } ;
                       % in future, add panel of new method here

%% function that create the method specific panel (Optional, for GUI)
% Specify the name of the function of the panel design 
function_EchoCombine_method_panel = {'sepia_handle_panel_EchoCombine_optimum_weights',...
                                     'sepia_handle_panel_EchoCombine_medi_nonlinear_fit',...
                                     };
                             % in future, add panel of new method here
                             
%% functions that generate and read BFR algorithm parameters (Optional, for GUI)
% Specify the name of the function that export/import parameters to/from
% output config file
config_EchoCombine_function     = {'get_set_echocombine_optimumweigts',...
                                   'get_set_echocombine_medinonlinearfit',...
                                   };
                       % in future, add panel of new method here
