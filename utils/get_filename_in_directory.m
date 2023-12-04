%% [file,numFiles] = get_filename_in_directory(inputDir,pattern,ext)
%
% Input
% --------------
% inputDir      : input directory
% pattern       : string pattern in the filename
% ext           : file extension
%
% Output
% --------------
% file          : structure contains all filenames with the specified
%                 pattern and extension
% numFiles      : number of files satisfied the criteria
%
% Description: get all the filenames in a directory with specific string
% pattern and extension
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created:  20 Jan 2021
% Date modified:
%
%
function [file,numFiles] = get_filename_in_directory(inputDir,pattern,ext)

% get all files with specified extension
fileList  = dir(fullfile(inputDir,['*' ext '*']));

% initiate output
file        = struct();
% initiate counter
numFiles    = 0;    

% check all files in the list
for klist = 1:length(fileList)                      % ignore hidden file
    if ContainName(fileList(klist).name,pattern) && ~strcmp(fileList(klist).name(1),'.')
        numFiles = numFiles + 1;
        file(numFiles).name = fullfile(fileList(klist).folder,fileList(klist).name);
    end
end

end