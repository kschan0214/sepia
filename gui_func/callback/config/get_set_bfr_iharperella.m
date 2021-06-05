%% function output = function_name(input)
%
% Usage:
%
% Input
% --------------
%
% Output
% --------------
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 
% Date last modified:
%
%
function get_set_bfr_iharperella(h,mode,input)


str_pattern = {'.bfr.iteration'};

action_handle = {h.bkgRemoval.iHARPERELLA.edit.maxIter};

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