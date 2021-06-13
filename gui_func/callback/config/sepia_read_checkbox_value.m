%% sepia_read_checkbox_value(config_txt, str_pattern, action_handle)
%
% Input
% --------------
% config_txt    : variable contains config text
% str_pattern   : string pattern to be printed after algorParam parameter
% action_handle : handle of the GUI popup
%
% Output
% --------------
% val           : (optional) value of the field
% 
% Description: read checkbox value from config file to GUI
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 12 JUne 2021 (v1.0)
% Date modified:
%
%
function val = sepia_read_checkbox_value(config_txt, str_pattern, action_handle)

% get value between = and ;, should be either 1 or 0
val = get_num_as_string(config_txt, str_pattern, '=', ';');

% set the value on the GUI
set_non_nan_value(action_handle, 'Value', str2double(val));

end