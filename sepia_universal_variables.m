%% sepia_universal_variables
%
% Description: a script contains all the method names and related function
% names
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 6 March 2020
% Date modified: 6 May 2021 (v0.8.1.1)
% Date modified: 7 June 2021 (v1.0)
% Date modified: 4 August 2021 (v1.0.1)
%
% DO NOT change the variable name
% DO NOT change the order of the entities, add a new one at the end instead
%
%% Version
SEPIA_version = 'v1.1dev';

%% PATH
SEPIA_HOME = fileparts(mfilename('fullpath'));

%% General parameterss
gyro = 42.57747892; % Larmor frequency of 1H, in MHz/T

%% Total field recovery related parameters
% Echo combination method available in SEPIA
sepia_configuration_EchoCombine;
% methodEchoCombineName   = {'Optimum weights',...
%                            'MEDI nonlinear fit',...
%                            'MEDI nonlinear fit (Bipolar, testing)'};

% Methods to exclude unreliable voxels                       
methodExcludedName      = {'Weighting map',...
                           'Brain mask'};

sepia_configuration_unwrap

%% background field removal related parameters
sepia_configuration_BFR

% Method to remove residual B1 field
methodRefineName        = {'3D Polynomial',...
                           'Spherical harmonic',...
                           'None'};

%% QSM related parameters
sepia_configuration_QSM

% Referencing region available in SEPIA
tissueName              = {'None',...
                           'Brain mask',...
                           'CSF'};

%% SWI/SMWI realted parameters
sepia_configuration_SWISMWI

%% add-ons capability
sepia_load_addons
