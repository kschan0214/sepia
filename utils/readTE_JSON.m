%% TE = readTE_JSON(filenames)
%
% Input
% --------------
% filenames     : json file fullname(s)
%
% Output
% --------------
% TE            : echo times
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 13 June 2018
% Date modified: 12 April 2019
%
%
function [TE] = readTE_JSON(filenames)

if ~iscell(filenames)
    filenames = {filenames};
end

if exist('jsondecode','builtin') == 5
    isJSONDecode = true;
    json_header = jsondecode(fileread(filenames{1}));
    conversion_software = json_header.ConversionSoftware;
else
    isJSONDecode = false;
    conversion_software = get_value_from_file(filenames{1},'ConversionSoftware');
end

TE = [];
if ~isempty(strfind(conversion_software,'dcm2niix'))
    % expect mutiple JSON files with dcm2niix
    if isJSONDecode
        for kf = 1:length(filenames)
            json_header = jsondecode(fileread(filenames{kf}));
            TE = [TE, json_header.EchoTime];
        end
    else
        for kf = 1:length(filenames)
            TE = [TE,str2double(get_value_from_file(filenames{kf},'EchoTime'))];
        end
    end
elseif ~isempty(strfind(conversion_software,'dicm2nii'))
    % expect only one JSON file with dicm2nii
    if isJSONDecode
        json_header = jsondecode(fileread(filenames{1}));
        TE = json_header.EchoTimes;
    else
        fid = fopen(filenames{1},'r');
        % read file line by line
        line = fgetl(fid);
        % start reading lines
        isNextline = false;
        while ischar(line)
            
            if isNextline
                % find chars ': '
                if ~isempty(strfind(line,']'))
                    separator = ']';
                    isNextline = false;
                else
                    separator = ',';
                end
                k_b = strfind(line,separator);
                % account for two chars ':' and ' '
                TE = [TE,str2double(line(1:k_b-1))];
            end

            % look for string 'Echo time'
            if ~isempty(strfind(line, 'EchoTimes'))
                isNextline = true;
            end
            
            % start the next line
            line = fgetl(fid);
        end
        fclose(fid);
    end
    
else
    error('Only dcm2niix and dicm2nii JSON files are supported.');

end

TE = unique(TE,'sorted');
% TE = sort(TE,'ascend');

end

%% Get value of the tag 
function str=get_str(list_info)
    % find chars ': '
    k_b = strfind(list_info,': ');
    % account for two chars ':' and ' '
    str=(list_info(k_b(1)+2:end));
end

function target_value = get_value_from_file(fname,target_string)

% read mode
fid = fopen(fname,'r');
% read file line by line
line = fgetl(fid);
    
% start reading lines
while ischar(line)
    % look for target string
    if ~isempty(strfind(line, target_string))
        target_value = get_str(line);
        break
    end
    % start the next line
    line = fgetl(fid);
end

fclose(fid);
end