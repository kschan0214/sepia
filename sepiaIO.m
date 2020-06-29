%% sepiaIO(input,output,maskFullName,algorParam)
%
% Input
% --------------
% input         :   input directory contains NIfTI files or structure containing filenames  
% output        :   output directory that stores the output
% maskFullName  :   mask filename
% algorParam    :   structure contains method and method specific parameters
%
% Description: This is a I/O level wrapper to coordinate SEPIA processing
% pipeline based on input algorithm parameter (sub)structure
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 29 June 2020 (v0.8.0)
% Date modified: 
%
%
function sepiaIO(input,output,maskFullName,algorParam)

% check if the input algorithm parameters contain any specific tasks
isUnwrapParam   =  isfield(algorParam,'unwrap');
isBRFParam      =  isfield(algorParam,'bfr');
isQSMParam      =  isfield(algorParam,'qsm');

% input algorithm parameters must have at least one specific task
isValid = or(or(isUnwrapParam,isBRFParam),isQSMParam);
if ~isValid
    error(['Input ''algorParam'' contains no parameter for any tasks in SEPIA.',...
           'Please check your pipeline configuration file.']);
end

% One-stop processing contains all task parameter
isOneStop = and(and(isUnwrapParam,isBRFParam),isQSMParam);

% determine which pipeline the data will go
if isOneStop
    % if algorParam contains all three task, then execute one-stop
    % processing pipeline
    disp('Running one-stop processing pipeline...');
    SepiaIOWrapper(input,output,maskFullName,algorParam);
else
    if and(isUnwrapParam,isBRFParam) || and(isUnwrapParam,isQSMParam) || and(isBRFParam,isQSMParam)
        % pipeline with 2 tasks is not allowed. Either 3 tasks or 1 task is
        % allowed
        error(['Input ''algorParam'' contains parameters for two tasks.',...
               'It can only contain parameters for either all three tasks (One-stop processing) or a single task.',...
               'Please check your pipeline configuration file.']);
    else
        % single task application pipeline
        if isUnwrapParam
            disp('Running total field recovery and phase unwrapping processing pipeline...');
            UnwrapPhaseMacroIOWrapper(input,output,maskFullName,algorParam);
        elseif isBRFParam
            disp('Running background field removal processing pipeline...');
            BackgroundRemovalMacroIOWrapper(input,output,maskFullName,algorParam);
        elseif isQSMParam
            disp('Running QSM dipole inversion processing pipeline...');
            QSMMacroIOWrapper(input,output,maskFullName,algorParam);
        end
    end
end