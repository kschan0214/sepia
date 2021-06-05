%% get_set_qsm_myQSM(h,mode,input)
%
% Input
% --------------
% h             : structure contains all handles of SEPIA
% mode          : 'set' - extract information from GUI to config file
%                 'get' - extract information from config file to GUI
% input         : if mode is 'set' then input should be fid
%                 if mode is 'get' then input should be confige file text
%
% Description: Information communication between config file and GUI
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 6 March 2020
% Date last modified:
%
%
function get_set_qsm_myQSM(h,mode,input)

str_pattern = {'.qsm.threshold'};

action_handle = {h.qsm.myQSM.edit.threshold};

switch lower(mode)
    case 'set'
        fid = input;
        
        fprintf(fid,'algorParam%s = %s ;\n'	,str_pattern{1},get(action_handle{1},	'String'));
        
    case 'get'
        
        config_txt = input;
        
        % first edit field
        val             = get_num_as_string(config_txt, str_pattern{1}, '=', ';');
        set_non_nan_value(action_handle{1},'String',val)

end