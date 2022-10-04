%% get_set_r2s_NLLS(h,mode,input)
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
function get_set_r2s_NLLS(h,mode,input)

str_pattern = {'.r2s.isParallel',...
               };

action_handle = {h.r2s.nlls.checkbox.isParallel,...
                 };
           
% menuFitType = {'Magnitude','Complex','Mixed'};

switch lower(mode)
    case 'set'
        fid = input;
        
        k = 1;
        fprintf(fid,'algorParam%s = %i ;\n'             ,str_pattern{k},get(action_handle{k},	'Value'));
        
    case 'get'
        
        config_txt = input;

        k = 1;
        pattern_curr    = str_pattern{k};
        val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
        set_non_nan_value(action_handle{k}, 'Value', str2double(val))

end