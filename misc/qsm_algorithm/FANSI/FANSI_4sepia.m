%% chi = FANSI_4sepia(phase,magn,spatial_res,alpha,mu,noise,options,B0_dir)
%
% Description: Wrapper to use FANSI under Sepia 
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 1 April 2019
% Date last modified:
%
%
function chi = FANSI_4sepia(phase,magn,spatial_res,alpha,noise,options,B0_dir)

params = [];

% additional parameters for sepia
params.maxOuterIter = options.iterations;
params.tol_update   = options.update;

% original FANSI function 
N = size(phase);
params.N = N;


if nargin > 5
    if isfield(options,'kernel_mode')
         kmode = options.kernel_mode;
    else
        kmode = 0;
    end
else
        kmode = 0;
end

% % JAC
% if nargin < 7
%     B0_dir = [0,0,1]; % Assuming slices perpendicular to B0 
% end

% params.K = dipole_kernel_fansi( N, spatial_res, kmode, B0_dir ); 
if nargin > 6
    params.K = dipole_kernel_angulated( N, spatial_res, B0_dir ); 
else
    params.K = dipole_kernel_fansi( N, spatial_res, kmode ); 
end

params.input = phase;
params.weight = magn; 
 
params.alpha1 = alpha;            % gradient L1 penalty
% params.mu1 = mu;                  % gradient consistency
params.mu1 = options.mu;

params.mu2 = options.mu2;


% if nargin > 6
%     if isfield(options,'gradient_mode')
%         gmode = options.gradient_mode;
%     else
%         gmode = 0;
%     end
% else
%         gmode = 0;
% end
% Gm = gradient_calc(magn,gmode); 
% Gm = max(Gm,noise); % Binary weighting not implemented here
% params.regweight = mean(Gm(:))./Gm;
if isfield(options,'gradient_mode')
    gmode = options.gradient_mode;
    Gm = gradient_calc(magn,gmode); 
    Gm = max(Gm,noise); % Binary weighting not implemented here
    params.regweight = mean(Gm(:))./Gm;
end



if nargin < 7
    options.nonlinear = true;
    options.tgv = false;
end

if options.isWeakHarmonic
    % weak harmonic regularisation
    params.beta = options.beta;
    params.muh  = options.muh;
    
    if options.nonlinear
        if options.tgv
            out = WH_nlTGV(params);
        else
            out = WH_nlTV_4sepia(params);
        end

    else
        if options.tgv
%             out = WH_wTGV_4sepia(params);
            out = WH_wTGV(params);
        else
            out = WH_wTV(params);
        end

    end
    
else
    
    if options.nonlinear
        if options.tgv
            out = nlTGV(params);
        else
            out = nlTV(params);
        end

    else
        if options.tgv
            out = wTGV(params);
        else
            out = wTV(params);
        end

    end
end

% Sepia check the size of the output
chi   = out.x;

end
