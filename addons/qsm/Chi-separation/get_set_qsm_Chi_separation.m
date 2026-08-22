%% get_set_qsm_Chi_separation(h,mode,input)
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
% Taechang Kim
% sakkar2@snu.ac.kr
% Date created: 01 December 2025
%
%
function get_set_qsm_Chi_separation(h,mode,input)

str_pattern = {'.qsm.solver',...
               '.qsm.R2s',...
               '.qsm.R2',...
               '.qsm.Dr',...
                };

action_handle = {h.qsm.Chi_separation.popup.solver,...
                 h.qsm.Chi_separation.edit.R2s,...
                 h.qsm.Chi_separation.edit.R2,...
                 h.qsm.Chi_separation.edit.Dr};
           
menuSolver       = {'Chi-separation-MEDI', 'Chi-separation-iLSQR', 'Chi-sepnet-R2*', 'Chi-sepnet-R2'''};

switch lower(mode)
    case 'set'
        fid = input;

        switch action_handle{1}.String{action_handle{1}.Value,1}
            case menuSolver{1}
                fprintf(fid,'algorParam%s = ''%s'' ;\n'     ,str_pattern{1},action_handle{1}.String{action_handle{1}.Value,1});
                fprintf(fid,'algorParam%s = ''%s'' ;\n'	,str_pattern{2},get(action_handle{2},	'String'));
                fprintf(fid,'algorParam%s = ''%s'' ;\n'	,str_pattern{3},get(action_handle{3},	'String'));
                fprintf(fid,'algorParam%s = %s ;\n'	,str_pattern{4},get(action_handle{4},	'String'));                
            case menuSolver{2}
                fprintf(fid,'algorParam%s = ''%s'' ;\n'     ,str_pattern{1},action_handle{1}.String{action_handle{1}.Value,1});
                fprintf(fid,'algorParam%s = ''%s'' ;\n'	,str_pattern{2},get(action_handle{2},	'String'));
                fprintf(fid,'algorParam%s = ''%s'' ;\n'	,str_pattern{3},get(action_handle{3},	'String'));
                fprintf(fid,'algorParam%s = %s ;\n'	,str_pattern{4},get(action_handle{4},	'String'));                
            case menuSolver{3}
                fprintf(fid,'algorParam%s = ''%s'' ;\n'     ,str_pattern{1},action_handle{1}.String{action_handle{1}.Value,1});
                fprintf(fid,'algorParam%s = ''%s'' ;\n'	,str_pattern{2},get(action_handle{2},	'String'));
                fprintf(fid,'algorParam%s = %s ;\n'	,str_pattern{4},get(action_handle{4},	'String'));                
            case menuSolver{4}
                fprintf(fid,'algorParam%s = ''%s'''' ;\n'     ,str_pattern{1},action_handle{1}.String{action_handle{1}.Value,1});
                fprintf(fid,'algorParam%s = ''%s'' ;\n'	,str_pattern{2},get(action_handle{2},	'String'));
                fprintf(fid,'algorParam%s = ''%s'' ;\n'	,str_pattern{3},get(action_handle{3},	'String'));
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

                for k = 2:4
                    pattern_curr    = str_pattern{k};
                    val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
                    set_non_nan_value(action_handle{k},'String',val)
                end
                
            case menuSolver{2}
                set_non_nan_value(action_handle{k},'Value',2)
                
                for k = 2:4
                    pattern_curr    = str_pattern{k};
                    val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
                    set_non_nan_value(action_handle{k},'String',val)
                end

            case menuSolver{3}
                set_non_nan_value(action_handle{k},'Value',3)
                
                for k = [2,4]
                    pattern_curr    = str_pattern{k};
                    val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
                    set_non_nan_value(action_handle{k},'String',val)
                end

            case menuSolver{4}
                set_non_nan_value(action_handle{k},'Value',4)
                
                for k = 2:4
                    pattern_curr    = str_pattern{k};
                    val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
                    set_non_nan_value(action_handle{k},'String',val)
                end
        end

end