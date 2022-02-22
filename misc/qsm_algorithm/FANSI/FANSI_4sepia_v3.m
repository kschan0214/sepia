function out = FANSI_4sepia_v3( phase, magn, alpha, options )
% FANSI wildcard function. This may be used to call the included functions
% globally, for a simplified use.
%
% Parameters: 
% phase: local field map data
% magn: magnitude data
% alpha: gradient L1 penalty, regularization weight
% options (structure with optional fields):
% options.isNonlinear: linear or nonlinear algorithm? (default = true, i.e nonlinear method)
% options.isTGV: TV or TGV regularization? (default = false, i.e. TV regularization)
% options.mu: ADMM Lagrange multiplier, regularization term.
% options.iterations: Maximum number of iterations
% options.update: Normalized update threshold to stop the algorithm.
% options.isGPU: Use GPU acceleration?
% options.gradientMode: this creates a gradient field to be used as
%                       regularization weights for TV/TGV, with:
%                       0 to use the vector field. 
%                       1 for the L1 norm, and 
%                       2 for the L2 norm
% options.noise: noise standard deviation in the complex signal, required
%                for the regularization weight based on the gradients.
% options.voxelSize: Spatial resolution vector, in mm, or normalized (mean = 1).
% options.B0_dir: main field direction, e.g. [0 0 1] (only for continuous
%                 kernel)
% options.kernelMode: dipole kernel formulation, with:
%                     0 for the continuous kernel proposed by Salomir, et al. 2003.
%                     1 for the discrete kernel proposed by Milovic, et al. 2017.
%                     2 for the Integrated Green function proposed by Jenkinson, et al. 2004
% Other function parameters are kept to their default values.
%
% Created by Carlos Milovic, 30.03.2017
% Modified by Julio Acosta-Cabronero, 26.05.2017
% Last modified by Carlos Milovic, 11.10.2021


% Parse arguments
params = [];

% Calculate the spatially variable regularization weight based on the
% gradient of the magnitude image (using R2* is also recommended).
if isfield(options,'gradientMode') % This variable was renamed in release 3.0
    gmode = options.gradientMode;
    Gm = gradient_calc(magn,gmode);
    if isfield(options,'noise')
        noise = options.noise;
        Gm = max(Gm,noise); % Binary weighting not implemented
    end
    params.regweight = mean(Gm(:))./Gm;
end


% Calcule the required dipole kernel
kmode = 0;
voxelSize = [1 1 1]; % This variable was renamed in release 3.0
if isfield(options,'voxelSize') % This variable was renamed in release 3.0
    voxelSize = options.voxelSize;
end

if isfield(options,'kernelMode') % This variable was renamed in release 3.0
    kmode = options.kernelMode;
end

if isfield(options,'B0_dir')
    params.K = dipole_kernel_angulated( size(phase), voxelSize, options.B0_dir );
else
    params.K = dipole_kernel_fansi( size(phase), voxelSize, kmode );
end


% Setup the required parameters
params.input = phase;
params.alpha1 = alpha;
params.weight = magn/max(magn(:));
% Optional parameters
if isfield(options,'mu')
    params.mu1 = options.mu;
end
if isfield(options,'isGPU') % This variable is new in release 3.0
    params.isGPU = options.isGPU;
end
if isfield(options,'iterations')
    params.maxOuterIter = options.iterations;
end
if isfield(options,'update')
    params.tolUpdate = options.update; % This variable was renamed in release 3.0
end

% Launch the desired solver
isNonlinear = true; % This variable was renamed in release 3.0
isTGV = false; % This variable was renamed in release 3.0
if isfield(options,'isNonlinear') % This variable was renamed in release 3.0
    isNonlinear = options.isNonlinear;
end
if isfield(options,'isTGV') % This variable was renamed in release 3.0
    isTGV = options.isTGV;
end

if options.isWeakHarmonic
    % weak harmonic regularisation
    params.beta = options.beta;
    params.muh  = options.muh;
    
    if isNonlinear
        if isTGV
            out = WH_nlTGV_4sepia_v3(params);
        else
            out = WH_nlTV(params);
        end

    else
        if isTGV
            out = WH_wTGV(params);
        else
            out = WH_wTV(params);
        end

    end
else
    
    if isNonlinear
        if isTGV
            out = nlTGV(params);
        else
            out = nlTV(params);
        end
    else
        if isTGV
            out = wTGV(params);
        else
            out = wTV(params);
        end
    end
end



end

