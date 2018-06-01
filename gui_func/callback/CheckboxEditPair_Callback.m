%% CheckboxEditPair_Callback(source,eventdata,handleToBeDisable,trueValue)
%
% Input
% --------------
% handleToBeDisable	: uicontrol to be disabled when the paired checkbox is
%                     chekced
% trueValue         : the condition of when the uicontrol should be
%                     disabled
%
% Description: This is a callback function to disable an edit field when a
%              checkbox is checked
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 1 June 2018
% Date last modified: 
%
%
function CheckboxEditPair_Callback(source,eventdata,handleToBeDisable,trueValue)
    % compare source value to trueValue
    % e.g. if an edit field needed to be disable by checking an checkbox
    % then trueValue of the checkbox is 0
    if source.Value == trueValue
        % if source is equal to trueValue then enables target handle
        set(handleToBeDisable,'Enable','on');
    else
        % if source do not equal to trueValue then disables target handle 
        set(handleToBeDisable,'Enable','off');
    end

end