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
function get_set_qsm_fansi(h,mode,input)

str_pattern = {'.qsm.tol',...
               '.qsm.maxiter',...
               '.qsm.lambda',...
               '.qsm.mu1',...
               '.qsm.mu2',...
               '.qsm.solver',...
               '.qsm.constraint',...
               '.qsm.gradient_mode',...
               '.qsm.isWeakHarmonic',...
               '.qsm.beta',...
               '.qsm.muh'};

action_handle = {h.qsm.FANSI.edit.tol,...
                 h.qsm.FANSI.edit.maxIter,...
                 h.qsm.FANSI.edit.lambda,...
                 h.qsm.FANSI.edit.mu,...
                 h.qsm.FANSI.edit.mu2,...
                 h.qsm.FANSI.popup.solver,...
                 h.qsm.FANSI.popup.constraints,...
                 h.qsm.FANSI.popup.gradientMode,...
                 h.qsm.FANSI.checkbox.isWeakHarmonic,...
                 h.qsm.FANSI.edit.beta,...
                 h.qsm.FANSI.edit.muh};
           

switch lower(mode)
    case 'set'
        fid = input;
        
        for k = 1:5
            fprintf(fid,'algorParam%s = %s ;\n'         ,str_pattern{k},get(action_handle{k},	'String'));
        end
        
        for k = 6:8
            fprintf(fid,'algorParam%s = ''%s'' ;\n'     ,str_pattern{k},action_handle{k}.String{action_handle{k}.Value,1});
        end
       
        k = k+1;
        fprintf(fid,'algorParam%s = %i ;\n'             ,str_pattern{k},get(action_handle{k},	'Value'));
        if get(action_handle{k},	'Value')
            for kk = k+1:k+2
                fprintf(fid,'algorParam%s = %s ;\n'    	,str_pattern{kk},get(action_handle{kk},	'String'));
            end
        end
        
    case 'get'
        
        config_txt = input;
        
        % first 5 edit fields
        for k = 1:5
            pattern_curr    = str_pattern{k};
            val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
            set_non_nan_value(action_handle{k},'String',val)
        end

        % solver
        k = k+1;
        pattern_curr    = str_pattern{k};
        val             = get_string_as_string(config_txt, pattern_curr);
        switch lower(val)
            case 'linear'
                set_non_nan_value(action_handle{k},'Value',2)
            case 'non-linear'
                set_non_nan_value(action_handle{k},'Value',1)
        end

        % constraint
        k = k+1;
        pattern_curr    = str_pattern{k};
        val             = get_string_as_string(config_txt, pattern_curr);
        switch lower(val)
            case 'tv'
                set_non_nan_value(action_handle{k},'Value',2)
            case 'tgv'
                set_non_nan_value(action_handle{k},'Value',1)
        end

        % gradient mode
        k = k+1;
        pattern_curr    = str_pattern{k};
        val             = get_string_as_string(config_txt, pattern_curr);
        switch lower(val)
            case 'vector field'
                set_non_nan_value(action_handle{k},'Value',1)
            case 'l1 norm'
                set_non_nan_value(action_handle{k},'Value',2)
            case 'l2 norm'
                set_non_nan_value(action_handle{k},'Value',3)
            case 'none'
                set_non_nan_value(action_handle{k},'Value',4)
        end

        % weak field harmonic
        k = k+1;
        pattern_curr    = str_pattern{k};
        val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
        set_non_nan_value(action_handle{k}, 'Value', str2double(val))
        % trigger popup callback to switch method panel
        feval(action_handle{k}.Callback{1},action_handle{k},[],{action_handle{k+1},action_handle{k+2}},1);

        for k = 10:11
            pattern_curr    = str_pattern{k};
            val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
            set_non_nan_value(action_handle{k},'String',val);
        end

end