%% EditInputMinMax_Callback(source,eventdata,isIntegerInput,lb,ub)
%
% Input
% --------------
% isIntegerInput	: true for returning integer, false for float value
% lb                : minimum value of the input
% ub                : maximum value of the input
%
% Description: This is a callback function to constraint the user input
%              value
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 1 June 2018
% Date last modified: 
%
%
function EditInputMinMax_Callback(source,eventdata,isIntegerInput,lb,ub)
% set the min/max input allowed in edit fields

    % check minimum
    if str2double(source.String)<lb
        source.String = num2str(lb);
    end
    
    % if maximum is set then check maximum
    if nargin==5
        if str2double(source.String)>ub
            source.String = num2str(ub);
        end
    end
    
    % make sure the input is interger for some fields
    if isIntegerInput
        source.String = num2str(round(str2double(source.String)));
    end
    
end