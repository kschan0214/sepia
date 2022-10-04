%% [r2s,t2s,m0] = R2sMacro(magnitude,TE,algorParam,mask,headerAndExtraData)
%
% Input
% --------------
% magnitude     : multi-echo data
% TE            : echo time
% algorParam    : structure contains fields with algorithm-specific parameter(s)
% mask          : (optional) signal mask
% headerAndExtraData : structure contains extra header info/data for the algorithm
%
% Output
% --------------
% r2s           : R2* map, in s^-1
% t2s           : T2* map, in s
% s0            : (T1w) signal intensity (a.u.)
%
% Description: This is a wrapper function to access individual dipole field
%              inversion algorithms for SEPIA (default: 'TKD')
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 4 October 2022
% Date modified: 
%
function [r2s,t2s,s0] = R2sMacro(magnitude,TE,algorParam,mask,headerAndExtraData)

if nargin<4
    mask = [];
end

sepia_universal_variables;
methodR2sName = lower(methodR2sName);

TE       = double(TE(:).');

algorParam          = check_and_set_SEPIA_algorithm_default(algorParam);
method              = algorParam.r2s.method;

headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);

if isempty(mask)
    dims = size(magnitude);
    mask = ones(dims(1:3));
end

disp('-----------');
disp('R2* mapping');
disp('-----------');
    
%% QSM algorithm
disp('Calculating T2*/R2*/S0 map...');
disp(['The following R2* algorithm will be used: ' method]);

% R2* mapping algorithm
for k = 1:length(wrapper_R2s_function)
    if strcmpi(method,methodR2sName{k})
        [r2s,t2s,s0] = feval(wrapper_R2s_function{k},magnitude,TE,mask,algorParam, headerAndExtraData);
    end
end


end

