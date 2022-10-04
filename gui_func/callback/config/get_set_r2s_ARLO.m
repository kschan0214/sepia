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
function get_set_r2s_ARLO(h,mode,input)

str_pattern = {'.r2s.s0mode',...
               };

action_handle = {h.r2s.arlo.popup.s0,...
                 };
           
menuS0 = {'1st echo','Weighted sum','Averaging'};

switch lower(mode)
    case 'set'
        fid = input;
        
        k = 1;
        fprintf(fid,'algorParam%s = ''%s'' ;\n'     ,str_pattern{k},action_handle{k}.String{action_handle{k}.Value,1});
        
    case 'get'
        
        config_txt = input;

        % solver
        k = 1;
        pattern_curr    = str_pattern{k};
        val             = get_string_as_string(config_txt, pattern_curr);
        switch lower(val)
            case lower(menuS0{1})
                set_non_nan_value(action_handle{k},'Value',1)
            case lower(menuS0{2})
                set_non_nan_value(action_handle{k},'Value',2)
            case lower(menuS0{3})
                set_non_nan_value(action_handle{k},'Value',3)
        end

end