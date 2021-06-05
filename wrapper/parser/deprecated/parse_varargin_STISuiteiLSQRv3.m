%% params = parse_varargin_STISuiteiLSQR(arg)
%
% Description: parser for QSM_iLSQR.p
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 6 September 2017
% Date last modified: 27 Feb 2020 (v0.8.0)
%
function params = parse_varargin_STISuiteiLSQRv3(arg)

% predefine parameters
params.niter=100;
params.tol_step1=0.01;
params.tol_step2=0.001;
params.Kthreshold=0.25;
params.padsize=[4,4,4];

% use user defined input if any
for kvar = 1:length(arg)
    if strcmpi(arg{kvar},'tol_step1')
        params.tol_step1 = double(arg{kvar+1});
        continue
    end
    if strcmpi(arg{kvar},'tol_step2')
        params.tol_step2 = double(arg{kvar+1});
        continue
    end
    if strcmpi(arg{kvar},'iteration')
        params.niter = double(arg{kvar+1});
        continue
    end
    if strcmpi(arg{kvar},'threshold')
        params.Kthreshold = double(arg{kvar+1});
        continue
    end
    if strcmpi(arg{kvar},'padsize')
        params.padsize = double(arg{kvar+1});
        continue
    end
end

end