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
% Date last modified: 1 March 2020 (v0.8.0)
%
%
function CheckboxEditPair_Callback(source,eventdata,handleToBeDisable,trueValue)

if ~iscell(handleToBeDisable)
    handleToBeDisable_cell{1} = handleToBeDisable;
else
    handleToBeDisable_cell = handleToBeDisable;
end
    % compare source value to trueValue
    % e.g. if an edit field needed to be disable by checking an checkbox
    % then trueValue of the checkbox is 0
    if source.Value == trueValue
        % if source is equal to trueValue then enables target handle
        for k = 1:length(handleToBeDisable_cell)
            set(handleToBeDisable_cell{k},'Enable','on');
        end
    else
        % if source do not equal to trueValue then disables target handle 
        for k = 1:length(handleToBeDisable_cell)
            set(handleToBeDisable_cell{k},'Enable','off');
        end
    end

end