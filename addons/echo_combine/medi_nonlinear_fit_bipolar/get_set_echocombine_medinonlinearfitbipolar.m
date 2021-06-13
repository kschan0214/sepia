%% get_set_echocombine_medinonlinearfitbipolar(h,mode,input)
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
% Date created: 12 June 2021 (v1.0)
% Date modified:
%
%
function get_set_echocombine_medinonlinearfitbipolar(h,mode,input)

sepia_universal_variables;

% these are the string to be printed in the pipeline config file
% for best practice, the entries should be matched with action_handle
str_pattern = {'.unwrap.unwrapMethod',...
               '',...
               '.unwrap.excludeMaskThreshold',...
               '.unwrap.excludeMethod'};

% these are the options available in the GUI
action_handle = {h.phaseUnwrap.MEDINonLinearfitBipolar.popup.phaseUnwrap,...
                 h.phaseUnwrap.MEDINonLinearfitBipolar.checkbox.excludeMask,...
                 h.phaseUnwrap.MEDINonLinearfitBipolar.edit.excludeMask,...
                 h.phaseUnwrap.MEDINonLinearfitBipolar.popup.excludeMethod};

switch lower(mode)
    case 'set'
        fid = input;
        
        % print spatial phase unwrapping method | popup 
        k = 1;
        sepia_print_popup_as_string(fid,str_pattern{k},action_handle{k});
        
        % more complicated operation for exclude mask options
        k = k+1;
        if get(action_handle{k},	'Value') % exclusion checkbox
            
            % exclusion threshold | edit
            k = k+1;
            sepia_print_edit_as_string(fid,str_pattern{k},action_handle{k});
            
            % exclusion method | popup
            k = k+1;
            sepia_print_popup_as_string(fid,str_pattern{k},action_handle{k});
        end
        

        
    case 'get'
        
        config_txt = input;
        
        % spatial phase unwrapping popup | popup
        k = 1;
        % change popup to the selected method
        sepia_read_popup_value(config_txt, str_pattern{k}, action_handle{k}, methodUnwrapName);
        % trigger popup callback to switch method panel
        % the input has to be matched with the callback function in the
        % sepia_handle_panel_EchoCombine_XXX.m
        feval(action_handle{k}.Callback{1},action_handle{k},[],h);
        
        % exclusion method | 'checkbox + edit + popup' combo
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
            
            % popup manu for thresholding method | popup
            k = k+1;
            sepia_read_popup_value(config_txt, str_pattern{k}, action_handle{k}, methodExcludedName);
            
        end

end