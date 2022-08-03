%% get_set_swismwi_swi_2dhamming(h,mode,input)
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
% Date created: 3 August 2022
% Date modified:
%
%
function get_set_swismwi_smwi(h,mode,input)

str_pattern = {'.swismwi.m',...
               '.swismwi.threshold',...
               '.swismwi.isParamagnetic',...
               '.swismwi.isDiamagnetic',...
               '.swismwi.ismIP',...
               '.swismwi.slice_mIP',...
               };

action_handle = {h.swismwi.smwi.edit.m,...
                 h.swismwi.smwi.edit.threshold,...
                 h.swismwi.smwi.checkbox.paramagnetic,...
                 h.swismwi.smwi.checkbox.diamagnetic,...
                 h.swismwi.smwi.checkbox.mIP,...
                 h.swismwi.smwi.edit.mIP,...
                 };

switch lower(mode)
    case 'set'
        fid = input;
        
        for k = 1:2
            fprintf(fid,'algorParam%s = %s ;\n'         ,str_pattern{k},get(action_handle{k},	'String'));
        end
        
        for k = 3:5
            fprintf(fid,'algorParam%s = %i ;\n'         ,str_pattern{k},get(action_handle{k},	'Value'));
        end
        if get(action_handle{k},	'Value')
            for kk = k+1
                fprintf(fid,'algorParam%s = %s ;\n'    	,str_pattern{kk},get(action_handle{kk},	'String'));
            end
        end
        
        
      % no load config operation supported at the moment
%     case 'get'
%         
%         config_txt = input;
%         
%         % Lambda
%         k = 1;
%         pattern_curr    = str_pattern{k};
%         val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
%         set_non_nan_value(action_handle{k},'String',val)
% 
%         % L-curve optimisation
%         k = k+1;
%         pattern_curr    = str_pattern{k};
%         val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
%         set_non_nan_value(action_handle{k}, 'Value', str2double(val))
%         % trigger popup callback to switch method panel
%         feval(action_handle{k}.Callback{1},action_handle{k},[],action_handle{k-1},0);

end