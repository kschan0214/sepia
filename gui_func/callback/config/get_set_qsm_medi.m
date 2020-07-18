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
function get_set_qsm_medi(h,mode,input)

str_pattern = {'.qsm.lambda',...
               '.qsm.wData',...
               '.qsm.zeropad',...
               '.qsm.percentage',...
               '.qsm.isSMV',...
               '.qsm.radius',...
               '.qsm.merit',...
               '.qsm.isLambdaCSF',...
               '.qsm.lambdaCSF'};

action_handle = {h.qsm.MEDI.edit.lambda,...
                 h.qsm.MEDI.edit.weightData,...
                 h.qsm.MEDI.edit.zeropad,...
                 h.qsm.MEDI.edit.percentage,...
                 h.qsm.MEDI.checkbox.smv,...
                 h.qsm.MEDI.edit.smv_radius,...
                 h.qsm.MEDI.checkbox.merit,...
                 h.qsm.MEDI.checkbox.lambda_csf,...
                 h.qsm.MEDI.edit.lambda_csf};
           
switch lower(mode)
    case 'set'
        fid = input;
        
        % lambda, wData, zeropad, percentage
        for k = 1:4
            fprintf(fid,'algorParam%s = %s ;\n'         ,str_pattern{k},get(action_handle{k},	'String'));
        end
        
        % isSMV
        k = k+1;
        fprintf(fid,'algorParam%s = %i ;\n'             ,str_pattern{k},get(action_handle{k},	'Value'));
        % radius
        k = k+1;
        fprintf(fid,'algorParam%s = %s ;\n'             ,str_pattern{k},get(action_handle{k},	'String'));
        % merit
        k = k+1;
        fprintf(fid,'algorParam%s = %i ;\n'             ,str_pattern{k},get(action_handle{k},	'Value'));
        % isLambdaCSF
        k = k+1;
        fprintf(fid,'algorParam%s = %i ;\n'             ,str_pattern{k},get(action_handle{k},	'Value'));
        % lambda_csf
        k = k+1;
        fprintf(fid,'algorParam%s = %s ;\n'             ,str_pattern{k},get(action_handle{k},	'String'));
        
    case 'get'
        
        config_txt = input;
        
        for k = 1:4
            pattern_curr    = str_pattern{k};
            val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
            set_non_nan_value(action_handle{k},'String',val)
        end

        % SMV
        k = k+1;
        pattern_curr    = str_pattern{k};
        val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
        set_non_nan_value(action_handle{k}, 'Value', str2double(val))
        % trigger popup callback to switch method panel
        feval(action_handle{k}.Callback{1},action_handle{k},[],action_handle{k+1},1);

        % SMV radius
        k = k+1;
        pattern_curr    = str_pattern{k};
        val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
        set_non_nan_value(action_handle{k},'String',val);


        % Merit
        k = k+1;
        pattern_curr    = str_pattern{k};
        val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
        set_non_nan_value(action_handle{k}, 'Value', str2double(val))

        % lambda CSF
        k = k+1;
        pattern_curr    = str_pattern{k};
        val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
        set_non_nan_value(action_handle{k}, 'Value', str2double(val))
        feval(action_handle{k}.Callback{1},action_handle{k},[],action_handle{k+1},1);

        % lambda CSF value
        k = k+1;
        pattern_curr    = str_pattern{k};
        val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
        set_non_nan_value(action_handle{k},'String',val);

end