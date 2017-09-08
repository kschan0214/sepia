% FANSI main function. This may be used to call the included functions
% globally, for a simplified use.
%
% Last modified by Carlos Milovic in 2017.03.30
%
function out = FANSI( phase, magn, spatial_res, alpha, mu, noise, options )
%
% input: 
% phase - local field map data
% magn - magnitude data
% spatial_res - Spatial resolution vector, in mm, or normalized (mean = 1).
% alpha - gradient L1 penalty, regularization weight
% mu - gradient consistency weight
% noise - noise standard deviation in the complex signal
% Optional:
% options.nonlinear - linear or nonlinear algorithm? (default = true)
% options.tgv - TV or TGV regularization? (default = false)
% options.kernel_mode - 0 for the continuous kernel proposed by Salomir, et al. 2003.
%                       1 for the discrete kernel proposed by Milovic, et al. 2017.
%                       2 for the Integrated Green function proposed by Jenkinson, et al. 2004
% options.gradient_mode - 0 to use the vector field. 
%                         1 for the L1 norm, and 
%                         2 for the L2 norm


params = [];
N = size(phase);
params.N = N;


if nargin > 6
    if isfield(options,'kernel_mode')
         kmode = options.kernel_mode;
    else
        kmode = 0;
    end
else
        kmode = 0;
end
params.K = dipole_kernel( N, spatial_res, kmode ); 

params.input = phase;
params.weight = magn; 
 
params.mu1 = mu;                  % gradient consistency
params.alpha1 = alpha;               % gradient L1 penalty


if nargin > 6
    if isfield(options,'gradient_mode')
         gmode = options.gradient_mode;
    else
        gmode = 0;
    end
else
        gmode = 0;
end
Gm = gradient_calc(magn,gmode); 
Gm = max(Gm,noise); % Binary weighting not implemented here
params.regweight = mean(Gm(:))./Gm;


if nargin < 7
    options.nonlinear = true;
    options.tgv = false;
end

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

