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
function get_set_qsm_iterTik(h,mode,input)

str_pattern = {'.qsm.solver',...
               '.qsm.threshold',...
               '.qsm.lambda',...
               '.qsm.tolerance'};

action_handle = {h.qsm.iterTik.popup.solver,...
                 h.qsm.iterTik.edit.threshold,...
                 h.qsm.iterTik.edit.lambda,...
                 h.qsm.iterTik.edit.tol};

menuSolver       = {'Truncated kspace division','Direct Tikhonov','Iterative Tikhonov'};

switch lower(mode)
    case 'set'
        fid = input;
        
        fprintf(fid,'algorParam%s = ''%s'' ;\n'     ,str_pattern{1},action_handle{1}.String{action_handle{1}.Value,1});
        
        switch action_handle{1}.String{action_handle{1}.Value,1}
            case menuSolver{1}
                fprintf(fid,'algorParam%s = %s ;\n'	,str_pattern{2},get(action_handle{2},	'String'));
            case menuSolver{2}
                fprintf(fid,'algorParam%s = %s ;\n'	,str_pattern{3},get(action_handle{3},	'String'));
            case menuSolver{3}
                fprintf(fid,'algorParam%s = %s ;\n'	,str_pattern{3},get(action_handle{3},	'String'));
                fprintf(fid,'algorParam%s = %s ;\n'	,str_pattern{4},get(action_handle{4},	'String'));
        end
        
    case 'get'
        
        config_txt = input;
        
        k = 1;
        pattern_curr    = str_pattern{k};
        val             = get_string_as_string(config_txt, pattern_curr);
        switch val
            case menuSolver{1}
                set_non_nan_value(action_handle{k},'Value',1)
                
                pattern_curr    = str_pattern{2};
                val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
                set_non_nan_value(action_handle{2},'String',val)
                
            case menuSolver{2}
                set_non_nan_value(action_handle{k},'Value',2)
                
                pattern_curr    = str_pattern{3};
                val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
                set_non_nan_value(action_handle{3},'String',val)
                
            case menuSolver{3}
                set_non_nan_value(action_handle{k},'Value',3)
                
                for k = 3:length(action_handle)
                    pattern_curr    = str_pattern{k};
                    val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
                    set_non_nan_value(action_handle{k},'String',val)
                end
        end
        feval(action_handle{1}.Callback{1},action_handle{1},[],h,menuSolver);

end