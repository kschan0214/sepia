%% sepia_read_popup_value(config_txt, str_pattern, action_handle, popup_list)
%
% Input
% --------------
% config_txt    : variable contains config text
% str_pattern   : string pattern to be printed after algorParam parameter
% action_handle : handle of the GUI popup
% popup_list    : popup list, in cell
%
% Description: read popup option from config file to GUI, without the need
%              to trigger any callback event
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 12 JUne 2021 (v1.0)
% Date modified:
%
%
function sepia_read_popup_value(config_txt, str_pattern, action_handle, popup_list)

% get option as string
val = get_string_as_string(config_txt, str_pattern);

if ~isnan(val)
    % matching popup list name
    for j = 1:length(popup_list)
        if strcmpi(val,popup_list{j})
            val = j;
            break
        end
    end

    % change popup manu
    set_non_nan_value(action_handle,'Value',val);
end

end