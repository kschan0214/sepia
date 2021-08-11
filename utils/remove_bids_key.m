%% str_new = remove_bids_key(str,key)
%
% Input
% --------------
% str           : BIDS compatible filename
% key           : key to be removed, string
%
% Output
% --------------
% str_new       : new filename without the specified key
%
% Description: remove a specific key from BIDS compatible filename
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 11 Augest 2021
% Date modified:
%
%
function str_new = remove_bids_key(str,key)

key = ['_' key '-'];

startIndex = strfind(str, key);

% key|value pair separator
startIndex_separator = strfind(str, '_');
if any(startIndex_separator > startIndex)
    % find the first separator after the key
    indAfter = startIndex_separator - startIndex;
    indAfter(indAfter<=0) = NaN;
    [~,indAfter] = min(abs(indAfter),[], 'omitnan');
else
    % in case the key located at the end of the filename
    startIndex_separator = strfind(str, '.');
    indAfter = startIndex_separator - startIndex;
    indAfter(indAfter<0) = NaN;
    [~,indAfter] = min(abs(indAfter),[], 'omitnan');
end

str_new = strcat(str(1:startIndex-1), str(startIndex_separator(indAfter):end));

end