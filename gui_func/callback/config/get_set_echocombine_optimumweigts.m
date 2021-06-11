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
function get_set_echocombine_optimumweigts(h,mode,input)

sepia_universal_variables;

% these are the string to be printed in the pipeline config file
str_pattern = {'.unwrap.unwrapMethod',...
               '.unwrap.isEddyCorrect',...
               '.unwrap.isSaveUnwrappedEcho',...
               '',...
               '.unwrap.excludeMaskThreshold',...
               '.unwrap.excludeMethod'};

% these are the options available in the GUI
action_handle = {h.phaseUnwrap.optimumWeights.popup.phaseUnwrap,...
                 h.phaseUnwrap.optimumWeights.checkbox.eddyCorrect,...
                 h.phaseUnwrap.optimumWeights.checkbox.saveEchoPhase,...
                 h.phaseUnwrap.optimumWeights.checkbox.excludeMask,...
                 h.phaseUnwrap.optimumWeights.edit.excludeMask,...
                 h.phaseUnwrap.optimumWeights.popup.excludeMethod};

switch lower(mode)
    case 'set'
        fid = input;
        
        % get popup selection
        k = 1;
        fprintf(fid,'algorParam%s = ''%s'' ;\n'     ,str_pattern{k},action_handle{k}.String{action_handle{k}.Value,1});
        
        % get checkbox value
        for k = 2:3
            fprintf(fid,'algorParam%s = %i ;\n'     ,str_pattern{k}, get(action_handle{k},'Value'));
        end
        
        % more complicated operation for exclude mask options
        k = k+1;
%         fprintf(fid,'algorParam%s = %i ;\n'             ,str_pattern{k},get(action_handle{k},	'Value'));
        if get(action_handle{k},	'Value')
            k = k+1;
            fprintf(fid,'algorParam%s = %s ;\n'    	,str_pattern{k},get(action_handle{k},	'String'));
            
            k = k+1;
            fprintf(fid,'algorParam%s = ''%s'' ;\n' ,str_pattern{k},action_handle{k}.String{action_handle{k}.Value,1});
        end
        

        
    case 'get'
        
        config_txt = input;
        
        % spatial phase unwrapping popup
        k = 1;
        pattern_curr    = str_pattern{k};
        val             = get_string_as_string(config_txt, pattern_curr);
        % method the user chosen will affect if exclusion method can be used or not 
        for j = 1:length(methodUnwrapName)
            if strcmpi(val,methodUnwrapName{j})
                % change method popup manu
                set_non_nan_value(action_handle{k},'Value',j)
                % trigger popup callback to switch method panel
                % this has to be matched with the callback function in the
                % sepia_handle_panel_EchoCombine_XXX.m
                feval(action_handle{k}.Callback{1},action_handle{k},[],h);
            end
        end
        
        % Bipolar correction and save unwrapped echo phase checkboxes
        for k = 2:3
            pattern_curr    = str_pattern{k};
            val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
            set_non_nan_value(action_handle{k}, 'Value', str2double(val))
        end
        
        % check for exclusion threshold
        k = k+2;    % there is no exclusion checkbox printed, so skip 1
        pattern_curr    = str_pattern{k};
        val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
        if isnan(val)
            % trigger checkbox
            set(action_handle{k-1},'Value',0);    % uncheck checkbox
            feval(action_handle{k-1}.Callback{1},action_handle{k-1},[],{action_handle{k},action_handle{k+1}},1);
        else
            % trigger checkbox
            set(action_handle{k-1},'Value',1);    % check checkbox
            feval(action_handle{k-1}.Callback{1},action_handle{k-1},[],{action_handle{k},action_handle{k+1}},1);
            
            % modifiy edit field value
            set_non_nan_value(action_handle{k}, 'String', val);
            
            % popup manu for thresholding method
            k = k+1;
            pattern_curr    = str_pattern{k};
            val      	= get_string_as_string(config_txt, pattern_curr);
            if ~isnan(val)
                % matching algorithm name
                for j = 1:length(methodExcludedName)
                    if strcmpi(val,methodExcludedName{j})
                        val = j;
                        break
                    end
                end

                % change method popup manu
                set_non_nan_value(action_handle{k},'Value',val);
            end
        end

end