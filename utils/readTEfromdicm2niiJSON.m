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
% Date created: 6 September 2018
% Date last modified:
%
%
% function [TE] = readTEfromJSON(filenames)
function [TE] = readTEfromdicm2niiJSON(filenames)

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
    
    % start reading lines
    while ischar(line)

        % look for string 'Echo time'
        if ~isempty(strfind(line, 'EchoTimes'))
            % in dcm2niix format echo time stored in unit of s, has to
            % convert to second
    %         TE(counter) = get_str(line);
            line = fgetl(fid);
            while ~isnan(str2double(line))
                TE = [TE,str2double(line)];
                line = fgetl(fid);
            end
            indicator_end = strfind(line, ']');
            TE = [TE,str2double(line(1:indicator_end-1))];
            break;
        end
        
        line = fgetl(fid);
            
    end

    fclose(fid);
end

TE = sort(TE,'ascend');
end
