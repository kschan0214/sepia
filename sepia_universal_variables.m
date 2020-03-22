%% sepia_universal_variables
%
% Description: a script contains all the method names and related function
% names
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 6 March 2020
% Date modified: 
%
% DO NOT change the variable name
% DO NOT change the order of the entities, add a new one at the end instead
%
%% General parameterss
gyro = 42.57747892; % Larmor frequency of 1H, in MHz/T

%% Total field recovery related parameters
% Echo combination method available in SEPIA
methodEchoCombineName   = {'Optimum weights',...
                           'MEDI nonlinear fit',...
                           'MEDI nonlinear fit (Bipolar, testing)'};

% Methods to exclude unreliable voxels                       
methodExcludedName      = {'Weighting map',...
                           'Brain mask'};

sepia_configuration_unwrap

%% background field removal related parameters
sepia_configuration_BFR

%% QSM related parameters
sepia_configuration_QSM

% Referencing region available in SEPIA
tissueName              = {'None',...
                           'Brain mask',...
                           'CSF'};
