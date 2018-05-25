%% TE = readTEfromText(filename)
%
% Input
% --------------
% filename      : text file fullname
%
% Output
% --------------
% TE            : echo times
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 22 April 2018
% Date last modified: 25 May 2018
%
%
function [TE,TR,FA,voxelSize] = readTEfromText(filename)
% read mode
fid = fopen(filename,'r');
% read file line by line
line = fgetl(fid);
% set counter of TE
counter = 1;

TE =[]; TR=[];

% start reading lines
while ischar(line)
    
    if ~isempty(strfind(line, 'Repetition time'))
        TR = get_str(line) * 1e-3;
    end
    
    if ~isempty(strfind(line, 'Flip angle'))
        FA = get_str(line) * 1e-3;
    end
    
    if ~isempty(strfind(line, 'Slice thickness'))
        vz = get_str(line);
    end
    
    if ~isempty(strfind(line, 'Voxel size x'))
        vx = get_str(line);
    end
    
    if ~isempty(strfind(line, 'Voxel size y'))
        vy = get_str(line);
    end
    
    % look for string 'Echo time'
    if ~isempty(strfind(line, 'Echo time'))
        % in MRIConvert format echo time stored in unit of ms, has to
        % convert to second
%         TE(counter) = get_str(line) * 1e-3;
        TE = [TE, get_str(line) * 1e-3];
        % prepare next TE
        counter = counter + 1;
    end
    % start the next line
    line = fgetl(fid);
end

voxelSize = [vx, vy, vz];

fclose(fid);

end

%% Get value of the tag 
function str=get_str(list_info)
    % find chars ': '
    k_b = strfind(list_info,': ');
    % account for two chars ':' and ' '
    str=str2double(list_info(k_b(1)+2:end));
end