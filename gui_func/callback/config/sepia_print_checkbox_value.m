%% sepia_print_checkbox_value(fid,str_pattern,action_handle)
%
% Input
% --------------
% fid           : fid of the sepia config file
% str_pattern   : string pattern to be printed after algorParam parameter
% action_handle : handle of the GUI popup
%
% Output
% --------------
% val           : (optional) output
%
% Description: print GUI checkbox value to sepia config file
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 12 June 2021 (v1.0)
% Date modified:
%
%
function val = sepia_print_checkbox_value(fid,str_pattern,action_handle)

val = get(action_handle,'Value');

fprintf(fid,'algorParam%s = %i ;\n'     ,str_pattern, val);

end