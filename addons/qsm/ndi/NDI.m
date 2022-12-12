%% chi = NDI(localField,mask,voxelSize,varargin)
%
% Usage: chi = NDI(localField,mask,voxelSize,'b0dir',B0_dir,'weight',weight,...
%                   'iteration',200,'stepsize',1,'tolerance',0.5);
%
% Input
% --------------
% localField    : local (tissue) field
% mask          : signal mask
% voxelSize     : spatial resolution (in mm)
% flag/value pairs:
% -----------------
% 'tol'            : tolerance RMSE
% 'iteration'      : maximum number of iterations
% 'stepsize'       : step size of gradient descent
% 'weight'         : weighting mask
% 'b0dir'          : direction of main field 
%
% Output
% --------------
%
% Description: Nonlinear Dipole inversion (NDI) 
%
% Modified based on script_ndi.m in https://github.com/polakd/NDI_Toolbox
%
% Referece: D Polak, I Chatnuntawech, J Yoon, S Srinivasan Iyer, J Lee, 
% K Setsompop, and B Bilgic. 
% VaNDI: Variational Nonlinear Dipole Inversion enables QSM without free 
% parameters (program number 0319). 
% Proceedings of the International Society for Magnetic Resonance in 
% Medicine 2019, Montreal Canada.
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 5 June 2019
% Date modified: 29 September 2022 (v1.1.1) add gpu compatibility
%
%
function chi = NDI(localField,mask,voxelSize,varargin)

% get matrix size
matrixSize = size(mask);

% check optional input and set default
[tolerance,stepSize,iteration,weight,b0dir,isGPU] = parse_varargin_NDI(varargin);

% verbose
disp('Non-linear Dipole Inversion (NDI)');
disp('---------------------------------');
disp('The following parameters are used:');
disp(['Tolerance: '             num2str(tolerance)]);
disp(['Max. itration: '         num2str(iteration)]);
disp(['Gradient step size: '    num2str(stepSize)]);

% if no weighting map than use signal mask
if isempty(weight)
    weight = mask;
end

% 20220930 KC: the difference is minor weight vs weight^2
weight = weight .^2;

% create dipole kernel
dipoleKernel = DipoleKernel(matrixSize,voxelSize,b0dir);

% initialise susceptibilit map and gradient
chi = zeros(matrixSize, 'like',localField);
grad_prev = zeros(matrixSize, 'like',localField);

if isGPU
try
    % add gpu compatibility
    weight          = gpuArray(weight);
    dipoleKernel    = gpuArray(dipoleKernel);
    chi             = gpuArray(chi);
    localField      = gpuArray(localField);
    grad_prev       = gpuArray(grad_prev);
    mask            = gpuArray(mask);
    isGPU           = true;
catch
    isGPU = false;
end
end

% case of one B0 direction
tic
for t = 1:iteration
    temp = weight .* sin(ifftn(dipoleKernel .* fftn(chi)) - localField);

    grad_f = 2 * sum(ifftn(dipoleKernel .* fftn(temp)), 4);

    chi = chi - stepSize * real(grad_f);

    update_grad = rmse(grad_prev(mask==1), grad_f(mask==1));
    
    if mod(t,5) == 0 || t == iteration || t == 1
        disp(['iter: ', num2str(t), '   grad update:', num2str(update_grad)])
    end

    if update_grad < tolerance
        if mod(t,5) ~= 0 || t ~= iteration || t ~= 1
            disp(['iter: ', num2str(t), '   grad update:', num2str(update_grad)])
        end
        break
    end

    grad_prev = grad_f;
end
toc

% masking output
chi = chi .* mask;

if isGPU
    chi = gather(chi);
end

end