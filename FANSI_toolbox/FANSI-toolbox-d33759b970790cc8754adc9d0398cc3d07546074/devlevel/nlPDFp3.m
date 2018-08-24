% QSM with a Total Variation regularization with spatially variable fidelity and regularization weights.
% This uses ADMM to solve the functional.
%
% Based on the code by Bilgic Berkin at http://martinos.org/~berkin/software.html
% Last modified by Carlos Milovic in 2017.03.30
%
function out = nlPDFp3(params)
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
    if isfield(params,'outermask')
         outermask = params.outermask;
    else
        outermask = mask;
    end
    

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
    
    if ~isfield(params,'delta_tol')
        delta_tol = 1e-6;
    else
        delta_tol = params.delta_tol;
    end
    
    
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


kernel = params.K;

if isfield(params,'precond')
    precond = params.precond;
else
    precond = true;
end

background = 0;
if precond
    
    
    
    if isfield(params,'ph_unwrap')
        z2  = params.ph_unwrap;
    else
        %[z2, lapunwrap, delta] = iterlap_unwrap( params.input, mask, weight, 250 );
        z2 = unwrap(params.input.*mask, [1 1 1]);
    end
    
    phase_unwrapped = z2;
    
    imagesc3d2(mask.*z2, N/2, 106, [90,90,-90], [-30,30], 0, 'z2input-iteruw');
    imagesc3d2(mask.*z2, N/2, 107, [90,90,-90], [-6,6], 0, 'z2-lineal');
    
    
    
    if isfield(params,'airmodel')
        phi_e  = params.airmodel;
    else
        phi_e = backmodel( outer_mask, 0.0, 2 );
    end

% LSqrt-Fit background model to unwrapped data, with linear functions
 [ background, kappa ] = fitmodels( z2, mask, weight, phi_e );
    
imagesc3d2(background, N/2, 100, [90,90,-90], [-6,6], 0, 'background');
        display(kappa);
        %display(phi0);
%     z2 = zeros(N,'single');    
% 
% z2 = params.input;
% 
% for i = 1:150
%     out_old = z2;
% z2 = z2 + 2*pi*round( (background - z2)/(2*pi) );
% 
% if sum(abs(out_old(:)-z2(:))) < 1
%     break;
% end
% 
% end
z2=phase_unwrapped-phi_e;%phase_unwrapped-background;
%[z2alt, z2, res] = iterlap_unwrap(angle(exp(1i*(params.input-background))), ones(N), weight, 250 );
%imagesc3d2(mask.*z2, N/2, 104, [90,90,-90], [-30,30], 0, 'z2');
z2 = mask.*z2;%mask.*unwrapLaplacian(z2,N,[1 1 1]);
%imagesc3d2(mask.*z2, N/2, 105, [90,90,-90], [-6,6], 0, 'z2');
%imagesc3d2(mask.*z2alt, N/2, 114, [90,90,-90], [-30,30], 0, 'z2a');
%imagesc3d2(mask.*z2alt, N/2, 115, [90,90,-90], [-6,6], 0, 'z2a');
        %display(max(z2(:)));
        %display(min(z2(:)));
        
        
    else
        z2 = zeros(N,'single');
    end
s2 = zeros(N,'single');
%return;
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
    
        
        
        dt = t/num_iter;
        
        rhs_z2 = mu2*real(ifftn(kernel.*fftn(z))+s2  );
        z2a =  rhs_z2 ./ (mu2) ;

        % Newton-Raphson method
       %if t > 100
        delta = inf;
        inn = 0;
        while (delta > delta_tol && inn < 50)
            inn = inn + 1;
            norm_old = norm(z2a(:));
            
            update = ( weight .* sin(z2a - params.input+background) + mu2*z2a - rhs_z2 ) ./ ( weight .* cos(z2a - params.input+background) + mu2 );            
        
            z2a = z2a - update;     
            delta = norm(update(:)) / norm_old;
        end        
        disp(delta)
       %else
           
           z2b =  (weight.*(phase_unwrapped-background)./(weight+mu2)) + rhs_z2./(weight + mu2);
           z2 = (1-dt)*z2b+dt*z2a;
       %end
        
            s2 = s2 + real(ifftn(kernel.*fftn(z))) - z2;
        
        
        s = s + x - z;  
    end
end
%x=x+kappa(5)*(1-outermask);
out.time = toc;toc

out.x = x;
out.phi = background+real(ifftn(kernel.*fftn((1-mask).*x)));


out.local = params.input;%(params.input - out.phi);
for i = 1:25
    out_old = out.local;
out.local = out.local + 2*pi*round( (out.phi - out.local)/(2*pi) );

if sum(abs(out_old(:)-out.local(:))) < 1
    break;
end

end
out.local = out.local-out.phi;
out.iter = t;



end
