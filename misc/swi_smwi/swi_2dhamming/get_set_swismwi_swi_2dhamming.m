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
function get_set_swismwi_swi_2dhamming(h,mode,input)

str_pattern = {'.swismwi.m',...
               '.swismwi.threshold',...
               '.swismwi.filterSize',...
               '.swismwi.isPositive',...
               '.swismwi.isNegative',...
               '.swismwi.ismIP',...
               '.swismwi.slice_mIP',...
               '.swismwi.echoCombineMethod',...
               };

action_handle = {h.swismwi.swi2dhamming.edit.m,...
                 h.swismwi.swi2dhamming.edit.threshold,...
                 h.swismwi.swi2dhamming.edit.filterSize,...
                 h.swismwi.swi2dhamming.checkbox.positive,...
                 h.swismwi.swi2dhamming.checkbox.negative,...
                 h.swismwi.swi2dhamming.checkbox.mIP,...
                 h.swismwi.swi2dhamming.edit.mIP,...
                 h.swismwi.swi2dhamming.popup.echoCombineMethod,...
                 };

switch lower(mode)
    case 'set'
        fid = input;
        
        for k = 1:3
            fprintf(fid,'algorParam%s = %s ;\n'         ,str_pattern{k},get(action_handle{k},	'String'));
        end
        
        for k = 4:6
            fprintf(fid,'algorParam%s = %i ;\n'         ,str_pattern{k},get(action_handle{k},	'Value'));
        end
        if get(action_handle{k},	'Value')
            for kk = k+1
                fprintf(fid,'algorParam%s = %s ;\n'    	,str_pattern{kk},get(action_handle{kk},	'String'));
            end
        end
        
        fprintf(fid,'algorParam%s = ''%s'' ;\n'     ,str_pattern{end},action_handle{end}.String{action_handle{end}.Value,1});
        
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