%% TE = readTEfromJSON(filenames)
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
% Date last modified:
%
%
function [TE] = readTEfromJSON(filenames)

if ~iscell(filenames)
    filenames = {filenames};
end

TE =[];
for kf = 1:length(filenames)
    filename = filenames{kf};
    
    % read mode
    fid = fopen(filename,'r');
    % read file line by line
    line = fgetl(fid);
    % set counter of TE
    counter = 1;
    
    % start reading lines
    while ischar(line)

        % look for string 'Echo time'
        if ~isempty(strfind(line, 'EchoTime'))
            % in dcm2niix format echo time stored in unit of s, has to
            % convert to second
    %         TE(counter) = get_str(line);
            TE = [TE, get_str(line)];
            % prepare next TE
            counter = counter + 1;
        end
        % start the next line
        line = fgetl(fid);
    end

    fclose(fid);
end

TE = sort(TE,'ascend');

end

%% Get value of the tag 
function str=get_str(list_info)
    % find chars ': '
    k_b = strfind(list_info,': ');
    % account for two chars ':' and ' '
    str=str2double(list_info(k_b(1)+2:end));
end