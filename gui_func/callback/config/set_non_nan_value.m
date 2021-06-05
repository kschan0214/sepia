%% set_non_nan_value(h,field,val)
%
% Input
% --------------
% h             : structure contains all SEPIA handles
% field         : handle field to be changed
% val           : new value
%
% Description: only set the parameter when it is valid
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 6 March 2020
% Date last modified:
%
% 
function set_non_nan_value(h,field,val)
    if ~isnan(val)
        set(h,field,val);
    end
end