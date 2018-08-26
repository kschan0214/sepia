% QSM with a Total Variation regularization with spatially variable fidelity and regularization weights.
% This uses ADMM to solve the functional.
%
% Based on the code by Bilgic Berkin at http://martinos.org/~berkin/software.html
% Last modified by Carlos Milovic in 2017.03.30
%
function out = nlPDFb(params)
%
% input: params - structure with the following required fields:
% params.input - local field map
% params.K - dipole kernel in the frequency space
% params.alpha1 - gradient L1 penalty, regularization weight
% params.mu1 - gradient consistency weight
%                 and the following optional fields:
% params.mu2 - fidelity consistency weight (recommended value = 1.0)
% params.maxOuterIter - maximum number of ADMM iterations (recommended = 50)
% params.tol_update - convergence limit, change rate in the solution (recommended = 1.0)
% params.weight - data fidelity spatially variable weight (recommended = magnitude_data). Not used if not specified
% params.regweight - regularization spatially variable weight. Not used if not specified
% params.N - array size
% params.precond - preconditionate solution (for stability)
%
% output: out - structure with the following fields:
% out.x - calculated susceptibility map
% out.iter - number of iterations needed
% out.time - total elapsed time (including pre-calculations)
%


tic

mu = params.mu1;
lambda = params.alpha1;
mask = params.mask;

    if isfield(params,'mu2')
         mu2 = params.mu2;
    else
        mu2 = 1.0;
    end
    
    if isfield(params,'N')
         N = params.N;
    else
        N = size(params.input);
    end

    if isfield(params,'maxOuterIter')
        num_iter = params.maxOuterIter;
    else
        num_iter = 50;
    end
    
    if isfield(params,'tol_update')
       tol_update  = params.tol_update;
    else
       tol_update = 1;
    end

    if isfield(params,'weight')
        weight = params.weight;
    else
        weight = ones(N);
    end
    weight = weight.*weight;
%     
%     if isfield(params,'regweight')
%         regweight = params.regweight;
%         if length(size(regweight)) == 3
%             regweight = repmat(regweight,[1,1,1,3]);
%         end
%     else
%         regweight = ones([N 3]);
%     end
    
    
    Wy = (weight.*params.input./(weight+mu2));
%     
% z_dx = zeros(N, 'single');
% z_dy = zeros(N, 'single');
% z_dz = zeros(N, 'single');
% 
% s_dx = zeros(N, 'single');
% s_dy = zeros(N, 'single');
% s_dz = zeros(N, 'single');

x = zeros(N, 'single');
z = zeros(N, 'single');
s = zeros(N, 'single');


    if isfield(params,'precond')
        precond = params.precond;
    else
        precond = true;
    end
    
    if precond
        z2 = params.input.*mask;%unwrapLaplacian(params.input.*mask);
    else
        z2 = zeros(N,'single');
    end
s2 = zeros(N,'single');

kernel = params.K;
% 
% 
% [k1, k2, k3] = ndgrid(0:N(1)-1,0:N(2)-1,0:N(3)-1);
% 
% E1 = 1 - exp(2i .* pi .* k1 / N(1));
% E2 = 1 - exp(2i .* pi .* k2 / N(2));
% E3 = 1 - exp(2i .* pi .* k3 / N(3));
% 
% E1t = conj(E1);
% E2t = conj(E2);
% E3t = conj(E3);
% 
% EE2 = E1t .* E1 + E2t .* E2 + E3t .* E3;
K2 = abs(kernel).^2;

%tic
for t = 1:num_iter
    % update x : susceptibility estimate
    
    x_prev = x;
    x = mu*(1-mask).*(z-s)./(lambda+mu*(1-mask));

    x_update = 100 * norm(x(:)-x_prev(:)) / norm(x(:));
    disp(['Iter: ', num2str(t), '   Update: ', num2str(x_update)])
    
    if x_update < tol_update
        break
    end
    
    if t < num_iter
        
        
        
            Fz = ( mu2*conj(kernel).*fftn(z2-s2) + mu*fftn((1-mask).*x+s) )./( mu2*K2+mu+eps );
            z = real(ifftn(Fz));
    
        
            z2 = Wy + mu2*real(ifftn(kernel.*fftn(z))+s2)./(weight + mu2);
        
            s2 = s2 + real(ifftn(kernel.*fftn(z))) - z2;
        
        
        s = s + x - z;  
    end
end
out.time = toc;toc

out.x = x;
out.phi = real(ifftn(kernel.*fftn((1-mask).*x)));
out.local = (params.input - out.phi);
out.iter = t;



end
