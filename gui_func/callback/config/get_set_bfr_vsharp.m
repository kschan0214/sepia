%% get_set_bfr_lbv(h,mode,input)
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
function get_set_bfr_vsharp(h,mode,input)


str_pattern = {'.bfr.radius'};

action_handle = {h.bkgRemoval.VSHARP.edit.maxRadius,...
                 h.bkgRemoval.VSHARP.edit.minRadius};

switch lower(mode)
    case 'set'
        fid = input;
        
        fprintf(fid,'algorParam%s = [%s:-1:%s] ;\n'	,str_pattern{1},get(action_handle{1},	'String'),get(action_handle{2},	'String'));
        
    case 'get'
        
        config_txt = input;
        
        % max radius
        pattern_curr    = str_pattern{1};
        val             = get_num_as_string(config_txt, pattern_curr, '[', ':');
        set_non_nan_value(action_handle{1},'String',val)

        % min radius
        val             = get_num_as_string(config_txt, pattern_curr, ':', ']');
        idx             = regexp(val,':');
        val             = val(idx+1:end);
        set_non_nan_value(action_handle{2},'String',val)

end