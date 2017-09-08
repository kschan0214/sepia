%% function output = function_name(input)
%
% Usage:
%
% Input
% --------------
%
% Output
% --------------
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 
% Date last modified:
%
%
function chi = qsmADMM(localField,mask,matrixSize,voxelSize,varargin)
% parse input argument
[mu,lambda,num_iter,tol_update] = parse_vararginADMM(varargin);
%% TV ADMM recon
z_dx = zeros(matrixSize, 'single');
z_dy = zeros(matrixSize, 'single');
z_dz = zeros(matrixSize, 'single');

s_dx = zeros(matrixSize, 'single');
s_dy = zeros(matrixSize, 'single');
s_dz = zeros(matrixSize, 'single');

x = zeros(matrixSize, 'single');

% dipole kernel
kernel = DipoleKernal(matrixSize,voxelSize);
K2 = abs(kernel).^2;

kspace = fftn(localField);
Dt_kspace = conj(kernel) .* kspace;

tic
for t = 1:num_iter
    % update x : susceptibility estimate
    tx = Ext .* fftn(z_dx - s_dx);
    ty = Eyt .* fftn(z_dy - s_dy);
    tz = Ezt .* fftn(z_dz - s_dz);
    
    x_prev = x;
    x = ifftn( (mu * (tx + ty + tz) + Dt_kspace) ./ (eps + K2 + mu * E2) );

    x_update = 100 * norm(x(:)-x_prev(:)) / norm(x(:));
    disp(['Iter: ', num2str(t), '   Update: ', num2str(x_update)])
    
    if x_update < tol_update
        break
    end
    
    if t < num_iter
        % update z : gradient varible
        Fx = fftn(x);
        x_dx = ifftn(Ex .* Fx);
        x_dy = ifftn(Ey .* Fx);
        x_dz = ifftn(Ez .* Fx);

        z_dx = max(abs(x_dx + s_dx) - lambda / mu, 0) .* sign(x_dx + s_dx);
        z_dy = max(abs(x_dy + s_dy) - lambda / mu, 0) .* sign(x_dy + s_dy);
        z_dz = max(abs(x_dz + s_dz) - lambda / mu, 0) .* sign(x_dz + s_dz);

        % update s : Lagrange multiplier
        s_dx = s_dx + x_dx - z_dx;
        s_dy = s_dy + x_dy - z_dy;            
        s_dz = s_dz + x_dz - z_dz;            
    end
end
toc

chi = real(x).*mask;


end

function [mu,lambda,num_iter,tol_update] = parse_vararginADMM(arg)
% predefine parameters
mu = 1e-2;              % gradient consistency
lambda = 2e-4;          % gradient L1 penalty
num_iter = 50;
tol_update = 1;

if ~isempty(arg)
    for kvar = 1:length(arg)
        if strcmpi(arg{kvar},'mu')
            mu = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'lambda')
            lambda = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'iteration')
            num_iter = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'tol')
            tol_update = arg{kvar+1};
        end
    end
end
end