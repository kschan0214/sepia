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
% Date last modified:
%
%
function TE = readTEfromText(filename)
% read mode
fid = fopen(filename,'r');
% read file line by line
line = fgetl(fid);
% set counter of TE
counter = 1;

TE =[];

% start reading lines
while ischar(line)
    % look for string 'Echo time'
    idx = strfind(line, 'Echo time');
    
    if ~isempty(idx)
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

fclose(fid);

end

%% Get value of the tag 
function str=get_str(list_info)
    % find chars ': '
    k_b = strfind(list_info,': ');
    % account for two chars ':' and ' '
    str=str2double(list_info(k_b(1)+2:end));
end