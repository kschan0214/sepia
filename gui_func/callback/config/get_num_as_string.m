%% str = get_num_as_string(A, str_pattern, start_indicator, end_indicator)
%
% Input
% --------------
% A             : confige file text
% str_pattern  	: string for matching in A
% 
% Output
% --------------
% str           : number as sting following str_pattern
%
% Description: extract parameter from config file to GUI
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 6 March 2020
% Date last modified:
%
%
function str = get_num_as_string(A, str_pattern, start_indicator, end_indicator)

% get the last position of thw string
str_end_idx = regexp(A,str_pattern,'end');

% check if there is more than one parameter with similar pattern
if length(str_end_idx) > 1
    for k = 1:length(str_end_idx)
        next_character = A(str_end_idx(k)+1);
        if isspace(next_character) || strcmp(next_character,'=') % if it is space to equal symbol
            str_end_idx = str_end_idx(k);
            break
        end
    end
end

if ~isempty(str_end_idx)

    indicatorS_idx  = regexp(A,start_indicator);
    indicatorE_idx 	= regexp(A,end_indicator);

    % get all characters between the indicators
    str = A(indicatorS_idx(find(indicatorS_idx > str_end_idx, 1 ))+1:indicatorE_idx(find(indicatorE_idx > str_end_idx, 1 ))-1);

    % remove all white space
    str = str(find(~isspace(str)));
else
    str = nan;
end

end
