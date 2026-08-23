%% isUpdated = update_toolbox_directory_entry(SEPIA_HOME, varName, newPath)
%
% Input
% --------------
% SEPIA_HOME    : SEPIA root directory (contains SpecifyToolboxesDirectory.m)
% varName       : name of the toolbox directory variable, e.g. 'FANSI_HOME'
% newPath       : full path to the toolbox directory
%
% Output
% --------------
% isUpdated     : true if SpecifyToolboxesDirectory.m was modified
%
% Description: comment out any existing (uncommented) assignment of
% varName in SpecifyToolboxesDirectory.m and append a new assignment
% pointing to newPath, but only if the current value actually differs
% from newPath (otherwise the file is left untouched).
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 23 August 2026
%
%
function isUpdated = update_toolbox_directory_entry(SEPIA_HOME, varName, newPath)

configFile = fullfile(SEPIA_HOME,'SpecifyToolboxesDirectory.m');

newPath = fullfile(fileparts(fullfile(newPath,filesep)),filesep);

% get the current value of varName, if any, from the config file
run(configFile);

if ~exist(varName,'var')
    isUpdateHome = true;
elseif isempty(eval(varName))
    isUpdateHome = true;
else
    curr_HOME    = fileparts(eval(varName));
    isUpdateHome = ~strcmp(curr_HOME, fileparts(newPath));
end

isUpdated = false;
if isUpdateHome

    fid             = fopen(configFile);
    directory_text  = textscan(fid, '%s', 'Delimiter','\n', 'CollectOutput',true);
    fclose(fid);
    lines = directory_text{1};

    % comment out any existing (uncommented) assignment of varName
    for j = 1:length(lines)
        if ContainName(lines{j}, lower(varName)) && ~strcmp(lines{j}(1),'%')
            lines{j} = ['% ' lines{j}];
        end
    end
    % append the new assignment
    lines{end+1} = sprintf('%s = ''%s'';', varName, newPath);

    fid = fopen(configFile, 'w');
    for j = 1:length(lines)
        fprintf(fid, '%s\n', lines{j});
    end
    fclose(fid);

    isUpdated = true;
end

end
