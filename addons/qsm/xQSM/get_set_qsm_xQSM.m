%% get_set_qsm_xQSM(h,mode,input)
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
% Date created: 24 September 2022
% Date last modified:
%
%
function get_set_qsm_xQSM(h,mode,input)

str_pattern = {'.qsm.solver',...
                };

action_handle = {h.qsm.xQSM.popup.solver,...
                };
           
menuSolver       = {'xQSM_invivo', 'xQSM_syn', 'xQSM_invivo_withNoiseLayer', 'Unet_invivo', 'Unet_syn'};

switch lower(mode)
    case 'set'
        fid = input;
        
        fprintf(fid,'algorParam%s = ''%s'' ;\n'     ,str_pattern{1},action_handle{1}.String{action_handle{1}.Value,1});
        
    case 'get'
        
        config_txt = input;
        
        k = 1;
        pattern_curr    = str_pattern{k};
        val             = get_string_as_string(config_txt, pattern_curr);
        switch val
            case menuSolver{1}
                set_non_nan_value(action_handle{k},'Value',1)
                
            case menuSolver{2}
                set_non_nan_value(action_handle{k},'Value',2)
                
            case menuSolver{3}
                set_non_nan_value(action_handle{k},'Value',3)
                
            case menuSolver{4}
                set_non_nan_value(action_handle{k},'Value',4)
                
            case menuSolver{5}
                set_non_nan_value(action_handle{k},'Value',5)
        end

end