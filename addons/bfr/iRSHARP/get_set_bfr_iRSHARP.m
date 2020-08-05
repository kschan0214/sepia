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
% Date created: 4 August 2020
% Date last modified:
%
%
function get_set_bfr_iRSHARP(h,mode,input)

str_pattern = {'.bfr.radius',...
               '.bfr.threshold',...
               '.bfr.constant'};

action_handle = {h.bkgRemoval.iRSHARP.edit.radius,...
                 h.bkgRemoval.iRSHARP.edit.threshold,...
                 h.bkgRemoval.iRSHARP.edit.constant};


switch lower(mode)
    case 'set'
        fid = input;
        
        for k = 1:length(action_handle)
            fprintf(fid,'algorParam%s = %s ;\n'	,str_pattern{k},get(action_handle{k},	'String'));
        end
        
    case 'get'
        
        config_txt = input;
        
        % first 3 edit fields
        for k = 1:length(action_handle)
            pattern_curr    = str_pattern{k};
            val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
            set_non_nan_value(action_handle{k},'String',val)
        end

end