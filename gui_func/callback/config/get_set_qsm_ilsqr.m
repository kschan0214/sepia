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
function get_set_qsm_ilsqr(h,mode,input)

str_pattern = {'.qsm.tol',...
               '.qsm.maxiter',...
               '.qsm.lambda',...
               '.qsm.optimise'};

action_handle = {h.qsm.iLSQR.edit.tol,...
                 h.qsm.iLSQR.edit.maxIter,...
                 h.qsm.iLSQR.edit.lambda,...
                 h.qsm.iLSQR.checkbox.lambda};
           

switch lower(mode)
    case 'set'
        fid = input;
        
        for k = 1:length(action_handle)-1
            fprintf(fid,'algorParam%s = %s ;\n'	,str_pattern{k},get(action_handle{k},	'String'));
        end
        k = k+1;
        fprintf(fid,'algorParam%s = %i ;\n'	,str_pattern{k},get(action_handle{k},	'Value'));
        
    case 'get'
        
        config_txt = input;
        
        for k = 1:length(action_handle)-1
            pattern_curr    = str_pattern{k};
            val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
            set_non_nan_value(action_handle{k},'String',val)
        end

        % L-curve optimisation
        k = k+1;
        pattern_curr    = str_pattern{k};
        val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
        set_non_nan_value(action_handle{k}, 'Value', str2double(val))
        % trigger popup callback to switch method panel
        feval(action_handle{k}.Callback{1},action_handle{k},[],action_handle{k-1},0);

end