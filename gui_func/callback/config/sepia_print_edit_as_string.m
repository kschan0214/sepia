%% sepia_print_edit_as_string(fid,str_pattern,action_handle)
%
% Input
% --------------
% fid           : fid of the sepia config file
% str_pattern   : string pattern to be printed after algorParam parameter
% action_handle : handle of the GUI popup
%
% Description: print GUI edit field value to sepia config file
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 12 June 2021 (v1.0)
% Date modified:
%
%
function sepia_print_edit_as_string(fid,str_pattern,action_handle)

fprintf(fid,'algorParam%s = %s ;\n'    	,str_pattern, get(action_handle,'String'));

end