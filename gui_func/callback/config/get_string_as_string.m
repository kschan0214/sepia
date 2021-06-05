%% str = get_string_as_string(A, str_pattern)
%
% Input
% --------------
% A             : confige file text
% str_pattern  	: string for matching in A
% 
% Output
% --------------
% str           : string following str_pattern
%
% Description: extract parameter from config file to GUI
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 6 March 2020
% Date last modified:
%
%
function str = get_string_as_string(A, str_pattern)

% get the last position of thw string
str_end_idx = regexp(A,str_pattern,'end');

if ~isempty(str_end_idx)

    indicator_idx  = regexp(A,'''');

    % get all characters between the indicators
    str = A(indicator_idx(find(indicator_idx > str_end_idx, 1 ))+1:indicator_idx(find(indicator_idx > str_end_idx, 1 )+1)-1);
else
    str = nan;
end

end