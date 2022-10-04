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
% Date modified: 27 Jan 2021 (v0.8.1)
%
%
function sepiaIO(input,output,maskFullName,algorParam)

%%% Step 1 %%%
currDir = pwd;
% 1.1: get and create output directory
output_index    = strfind(output, filesep);
outputDir       = output(1:output_index(end));
% if the output directory does not exist then create the directory
if exist(outputDir,'dir') ~= 7
    mkdir(outputDir);
end
% move to output directory
cd(outputDir)

% 1.2 log command window display to a text file
% use current time as unique identifier
identifier = datestr(datetime('now'),'yymmddHHMMSSFFF');

logFilename = fullfile(outputDir, ['run_sepia.log' identifier]);
while exist(logFilename,'file') == 2
    % update current time as unique identifier
    identifier = datestr(datetime('now'),'yymmddHHMMSSFFF');
    logFilename = fullfile(outputDir, ['run_sepia.log' identifier]);
end
diary(logFilename)

% check if any sepia config file exists in the output directory, if not
% then create one
check_and_create_sepia_config(input,output,maskFullName,algorParam,identifier);

% display the parent script
fn = dbstack('-completenames');
if length(fn) >=2 
    fprintf('Running script: %s\n',fn(2).file);
end

%%% Step 2 %%%
try 
% check if the input algorithm parameters contain any specific tasks
isUnwrapParam   =  isfield(algorParam,'unwrap');
isBRFParam      =  isfield(algorParam,'bfr');
isQSMParam      =  isfield(algorParam,'qsm');
isR2sParam      =  isfield(algorParam,'r2s');

% input algorithm parameters must have at least one specific task
isValid = or(or(or(isUnwrapParam,isBRFParam),isQSMParam),isR2sParam);
if ~isValid
    error(['Input ''algorParam'' contains no parameter for any tasks in SEPIA.',...
           'Please check your pipeline configuration file.']);
end

% One-stop processing contains all task parameter
isOneStop = and(and(isUnwrapParam,isBRFParam),isQSMParam);

%%% Step 3 %%%
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
        elseif isR2sParam
            disp('Running R2* mapping processing pipeline...');
            R2sIOWrapper(input,output,maskFullName,algorParam);
        end
    end
end

% turn off the log
diary off
% move back to original directory
cd(currDir)
    
catch ME
    
    % close log file
    disp('There was an error! Please check the command window/error message file for more information.');
    diary off
    % move back to original directory
    cd(currDir)
    
    % open a new text file for error message
    errorMessageFilename = fullfile(outputDir, ['run_sepia.error' identifier]);
    fid = fopen(errorMessageFilename,'w');
    fprintf(fid,'The identifier was:\n%s\n\n',ME.identifier);
    fprintf(fid,'The message was:\n\n');
    msgString = getReport(ME,'extended','hyperlinks','off');
    fprintf(fid,'%s',msgString);
    fclose(fid);
    
    % rethrow the error message to command window
    rethrow(ME);
end

%%
function check_and_create_sepia_config(input,output,maskFullName,algorParam,identifier)

sepia_universal_variables;

isConfigFileExist = false;

%%%%%%%%%%%% Step 1: Check if config file exists in the output dir %%%%%%%%%%%% 
output_index    = strfind(output, filesep);
outputDir       = output(1:output_index(end));

% create a new m file
configFileList    = dir(fullfile(outputDir,'*sepia_config*.m*'));
if ~isempty(configFileList)
    isConfigFileExist = true;
    disp('SEPIA configuration file already exists in the output directory.')
else   
    configFilename = fullfile(outputDir, ['sepia_config' identifier '.m']);
end

%%%%%%%%%%%% Step 2: Write config file %%%%%%%%%%%% 
if ~isConfigFileExist

fid = fopen(configFilename,'w');
% specify config file version
fprintf(fid,'%% This file is generated by SEPIA version: %s\n',SEPIA_version);
% general path
fprintf(fid,'%% add general Path\n');
fprintf(fid,'sepia_addpath\n\n');

fprintf(fid,'%% Input/Output filenames\n');

% 2.1: input data
if isstruct(input)
    fprintf(fid,'input(1).name = ''%s'' ;\n',input(1).name);
    fprintf(fid,'input(2).name = ''%s'' ;\n',input(2).name);
    fprintf(fid,'input(3).name = ''%s'' ;\n',input(3).name);
    fprintf(fid,'input(4).name = ''%s'' ;\n',input(4).name);
else
    fprintf(fid,'input = ''%s'' ;\n',input);
end

% 2.2: output
fprintf(fid,'output_basename = ''%s'' ;\n',output);

% 2.3: mask
fprintf(fid,'mask_filename = [''%s''] ;\n\n',maskFullName);

% 2.4: algorithm parameters
task_field = fieldnames(algorParam);
fprintf(fid,'algorParam=struct();\n');
% loop all the tasks
for k = 1:length(task_field)
    % task separator
    switch task_field{k}
        case 'general'
            fprintf(fid,'%% General algorithm parameters\n');
        case 'unwrap'
            fprintf(fid,'%% Total field recovery algorithm parameters\n');
        case 'bfr'
            fprintf(fid,'%% Background field removal algorithm parameters\n');
        case 'qsm'
            fprintf(fid,'%% QSM algorithm parameters\n');
    end
        
    curr_task_field = algorParam.(task_field{k});
    
    method_field = fieldnames(curr_task_field);
    % loop all method parameters
    for kk = 1:length(method_field)
        curr_method_value = curr_task_field.(method_field{kk});
        
        if ischar(curr_method_value)
            curr_method_value = sprintf('''%s''',curr_method_value);
        elseif ~isscalar(curr_method_value)
            curr_method_value = sprintf('[%s]',num2str(curr_method_value));
        elseif isnumeric(curr_method_value)
            curr_method_value = num2str(curr_method_value);
        end
        
        % export to config file
        fprintf(fid,'%s.%s.%s = %s;\n','algorParam',task_field{k},method_field{kk},curr_method_value);
    end
end

fprintf(fid,'\nsepiaIO(input,output_basename,mask_filename,algorParam);\n');

disp('No SEPIA configuration file exists in the output directory. The file is now generated.')
end