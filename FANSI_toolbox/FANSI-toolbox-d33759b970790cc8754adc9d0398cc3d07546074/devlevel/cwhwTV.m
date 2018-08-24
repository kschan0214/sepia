% Single-Step QSM with a Total Variation regularization with spatially variable fidelity and regularization weights.
% This uses ADMM to solve the functional.
%
% Last modified by Carlos Milovic in 2017.06.09
%
function out = cwhwTV(params)
%
% input: params - structure with the following required fields:
% params.input - local field map
% params.K - dipole kernel in the frequency space
% params.alpha1 - gradient L1 penalty, regularization weight
% params.mu1 - gradient consistency weight
%                 and the following optional fields:
% params.mu2 - fidelity consistency weight (recommended value = 1.0)
% params.beta - harmonic constrain weight
% params.muh - harmonic consistency weight (recommended value = beta)
% params.mask - ROI to calculate susceptibility values (if not provided, will be calculated from 'weight')
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
% out.phi - harmonic phase in radians
%


tic

mu1 = params.mu1;
alpha1 = params.alpha1;

    if isfield(params,'mu2')
         mu2 = params.mu2;
    else
        mu2 = 1.0;
    end
    
    if isfield(params,'mu')
         mu = params.mu;
    else
        mu = 1.0;
    end
    
    if isfield(params,'muh')
         muh = params.muh;
    else
        muh = 1.0;
    end
    if isfield(params,'beta')
         beta = params.beta;
    else
        beta = 1.0;
    end
    
    if isfield(params,'gamma')
         gamma = params.gamma;
    else
        gamma = 1.0;
    end
    
    if isfield(params,'mask')
         mask = params.mask;
    else
        mask = params.weight > 0;
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
    
    if isfield(params,'regweight')
        regweight = params.regweight;
        if length(size(regweight)) == 3
            regweight = repmat(regweight,[1,1,1,3]);
        end
    else
        regweight = ones([N 3]);
    end
    
    
    Wy = (weight.*params.input./(weight+mu2));
    
z1_dx = zeros(N, 'single');
z1_dy = zeros(N, 'single');
z1_dz = zeros(N, 'single');

s1_dx = zeros(N, 'single');
s1_dy = zeros(N, 'single');
s1_dz = zeros(N, 'single');

x = zeros(N, 'single');

phi_h = zeros(N, 'single');

z_h = zeros(N, 'single');
s_h = zeros(N, 'single');
z = zeros(N, 'single');
s = zeros(N, 'single');

    if isfield(params,'precond')
        precond = params.precond;
    else
        precond = true;
    end
    
    if precond
        z2 = Wy;
    else
        z2 = zeros(N,'single');
    end
s2 = zeros(N,'single');

kernel = params.K;


[k1, k2, k3] = ndgrid(0:N(1)-1,0:N(2)-1,0:N(3)-1);

E1 = 1 - exp(2i .* pi .* k1 / N(1));
E2 = 1 - exp(2i .* pi .* k2 / N(2));
E3 = 1 - exp(2i .* pi .* k3 / N(3));

E1t = conj(E1);
E2t = conj(E2);
E3t = conj(E3);

EE2 = E1t .* E1 + E2t .* E2 + E3t .* E3;
Lap = E1+E1t+E2+E2t+E3+E3t;
K2 = abs(kernel).^2;

eta1 = 10*mu2;
eta2 = 10*gamma;
% 
% [ky,kx,kz] = meshgrid((-N(2)/2):(N(2)/2-1), (-N(1)/2):(N(1)/2-1), (-N(3)/2):(N(3)/2-1));
% G = (kx).^2+(ky).^2+(kz).^2;
% G = fftshift( exp( -G/(10^2) ) );
% G = G/sum(G(:));
% phi_h = real(ifftn(fftn(G).*fftn(params.input.*mask)));
        
        r1 = zeros( N ,'single');
        r2 = zeros( N ,'single');
        y1 = zeros( N ,'single');
        y2 = zeros( N ,'single');
        z = zeros( N ,'single');
%tic
for t = 1:num_iter
    % update x : susceptibility estimate
    tx = E1t .* fftn(z1_dx - s1_dx);
    ty = E2t .* fftn(z1_dy - s1_dy);
    tz = E3t .* fftn(z1_dz - s1_dz);
    
    x_prev = x;
    %Dt_kspace = conj(kernel) .* fftn(z2-s2-phi_h);
    x = real(ifftn( (mu1 * (tx + ty + tz) + mu*fftn(z-s)) ./ (eps + mu + mu1 * EE2) ));

    x_update = 100 * norm(x(:)-x_prev(:)) / norm(x(:));
    disp(['Iter: ', num2str(t), '   Update: ', num2str(x_update)])
    
    if x_update < tol_update
        break
    end
    
        ll = alpha1/mu1;
    if t < num_iter
        % update z1 : gradient variable
        Fx = fftn(x);
        x_dx = real(ifftn(E1 .* Fx));
        x_dy = real(ifftn(E2 .* Fx));
        x_dz = real(ifftn(E3 .* Fx));
        
        z1_dx = max(abs(x_dx + s1_dx) - mask.*regweight(:,:,:,1)*ll, 0) .* sign(x_dx + s1_dx);
        z1_dy = max(abs(x_dy + s1_dy) - mask.*regweight(:,:,:,2)*ll, 0) .* sign(x_dy + s1_dy);
        z1_dz = max(abs(x_dz + s1_dz) - mask.*regweight(:,:,:,3)*ll, 0) .* sign(x_dz + s1_dz);
    
        % update s1 : Lagrange multiplier
        s1_dx = s1_dx + x_dx - z1_dx;
        s1_dy = s1_dy + x_dy - z1_dy;            
        s1_dz = s1_dz + x_dz - z1_dz;  
        
        
        
        z2 = Wy + mu2*real(ifftn(kernel.*fftn(mask.*z))+s2+phi_h )./(weight + mu2);
               
        
        phi_h = real(ifftn((muh * conj(Lap).*fftn(z_h-s_h) + mu2*fftn(z2-s2) + kernel.*fftn( gamma*(1-mask).*z-mu2*mask.*z ) ) ./ (eps + mu2 + muh * Lap.*conj(Lap)+gamma) ));
        
        z_h = muh*(real(ifftn(Lap.*fftn(phi_h)))+s_h)./(muh+beta*mask);
%         
%         r1 = zeros( N ,'single');
%         r2 = zeros( N ,'single');
%         y1 = zeros( N ,'single');
%         y2 = zeros( N ,'single');
%         z = zeros( N ,'single');
         % z subproblem - inner loop
%           for tz = 1:(num_iter/5)
%             z_prev = z;
            z = ( mu*(x+s) + eta1*mask.*(y1-r1)+eta2*(1-mask).*(y2-r2) )./( mu+eta1*mask+eta2*(1-mask) );
            
          
%           
%             z_update = 100 * norm(z(:)-z_prev(:)) / norm(z(:));
%     
%             if z_update < (2*tol_update)
%                break
%             end
%     
          
            y1 = real(ifftn( (mu2*conj(kernel).* fftn(z2-s2-phi_h)+eta1*fftn(mask.*z+r1))./(mu2*K2+eta1)  ));
            y2 = real(ifftn( (gamma*conj(kernel).* fftn(phi_h)+eta2*fftn((1-mask).*z+r2))./(gamma*K2+eta2)  ));
            
            r1 = r1 + mask.*z-y1;
            r2 = r2 + (1-mask).*z-y2;
%           end
        
        
            s = s + x - z;
        
        
        s2 = s2 + real(ifftn(kernel.*fftn(mask.*z))) - z2+phi_h;
        s_h = s_h + real(ifftn(Lap.*fftn(phi_h))) - z_h;
    end
end
out.time = toc;toc

out.x = x;
out.phi = phi_h;
out.iter = t;



end
