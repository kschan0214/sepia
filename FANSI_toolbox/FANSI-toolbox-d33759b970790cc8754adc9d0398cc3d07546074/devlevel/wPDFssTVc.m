% QSM with a Total Variation regularization with spatially variable fidelity and regularization weights.
% This uses ADMM to solve the functional.
%
% Based on the code by Bilgic Berkin at http://martinos.org/~berkin/software.html
% Last modified by Carlos Milovic in 2017.03.30
%
function out = wPDFssTVc(params)
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

mu1 = params.mu1;
alpha1 = params.alpha1;
mask = params.mask;
beta = params.beta;

    if isfield(params,'mu_e')
         mu_e = params.mu_e;
    else
        mu_e = beta;
    end
    
    if isfield(params,'mu_l')
         mu_l = params.mu_l;
    else
        mu_l = 1.0;
    end
    
    if isfield(params,'nue')
         nu_e = params.nu_e;
    else
        nu_e = 1.0;
    end
    
    if isfield(params,'nu_l')
         nu_l = params.nu_l;
    else
        nu_l = 1.0;
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
        regweight = mask;
        regweight = repmat(regweight,[1,1,1,3]);
        %regweight = ones([N 3]);
    end
    
    
    Wy = (weight.*params.input./(weight+mu_l));
    
z1_dx = zeros(N, 'single');
z1_dy = zeros(N, 'single');
z1_dz = zeros(N, 'single');

s1_dx = zeros(N, 'single');
s1_dy = zeros(N, 'single');
s1_dz = zeros(N, 'single');


x = zeros(N, 'single');
y_e = zeros(N, 'single');
r_e = zeros(N, 'single');
y_l = zeros(N, 'single');
r_l = zeros(N, 'single');


[ky,kx,kz] = meshgrid((-N(2)/2):(N(2)/2-1), (-N(1)/2):(N(1)/2-1), (-N(3)/2):(N(3)/2-1));
G = (kx).^2+(ky).^2+(kz).^2;
G = fftshift( exp( -G/(5^2) ) )./sum(G(:));
G = fftn(G);

    if isfield(params,'precond')
        precond = params.precond;
    else
        precond = true;
    end
    %params.input.*weight;
    if precond
        z_e = params.input.*weight./(weight+mu_e);%real(ifftn(G.*fftn(params.input.*params.weight)));
        z_l = zeros(N,'single');%(params.input-z_e).*weight./(weight+mu_l);
    else
        z_e = zeros(N,'single');
        z_l = zeros(N,'single');
    end
s_e = zeros(N,'single');
s_l = zeros(N,'single');

kernel = params.K;


[k1, k2, k3] = ndgrid(0:N(1)-1,0:N(2)-1,0:N(3)-1);

E1 = 1 - exp(2i .* pi .* k1 / N(1));
E2 = 1 - exp(2i .* pi .* k2 / N(2));
E3 = 1 - exp(2i .* pi .* k3 / N(3));

E1t = conj(E1);
E2t = conj(E2);
E3t = conj(E3);

EE2 = E1t .* E1 + E2t .* E2 + E3t .* E3;
K2 = abs(kernel).^2;

 a1 = weight+mu_l;
% a2 = weight;
 a3 = weight+mu_e/beta;
% detA = a1.*a3-a2.*a2+eps;
maskc = 1-mask;


%tic
for t = 1:num_iter
    % update x : susceptibility estimate
    
    x_prev = x;
    x = ( nu_l*mask.*(y_l-r_l) + nu_e*maskc.*(y_e-r_e) )./(nu_l*mask+nu_e*maskc+eps);
    
    x_update = 100 * norm(x(:)-x_prev(:)) / norm(x(:));
    disp(['Iter: ', num2str(t), '   Update: ', num2str(x_update)])
    
    if x_update < tol_update
        break
    end
    
        ll = alpha1/mu1;
    if t < num_iter
        
        y_e = real(ifftn( (mu_e*conj(kernel).*fftn(z_e-s_e) + nu_e*fftn(maskc.*x+r_e) )./( mu_e*K2+nu_e+eps ) ));
        
        tx = E1t .* fftn(z1_dx - s1_dx);
        ty = E2t .* fftn(z1_dy - s1_dy);
        tz = E3t .* fftn(z1_dz - s1_dz);
        y_l = real(ifftn( (mu_l*conj(kernel).*fftn(z_l-s_l) + nu_l*fftn(mask.*x+r_l) + mu1*(tx + ty + tz) )./( mu_l*K2+nu_l+mu1*EE2+eps ) ));
        
        r_e = r_e - y_e + maskc.*x;
        r_l = r_l - y_l + mask.*x;
        
        
        
        Fx = fftn(y_l);
        x_dx = real(ifftn(E1 .* Fx));
        x_dy = real(ifftn(E2 .* Fx));
        x_dz = real(ifftn(E3 .* Fx));
        
        z1_dx = max(abs(x_dx + s1_dx) - regweight(:,:,:,1)*ll, 0) .* sign(x_dx + s1_dx);
        z1_dy = max(abs(x_dy + s1_dy) - regweight(:,:,:,2)*ll, 0) .* sign(x_dy + s1_dy);
        z1_dz = max(abs(x_dz + s1_dz) - regweight(:,:,:,3)*ll, 0) .* sign(x_dz + s1_dz);
    
        % update s : Lagrange multiplier
        s1_dx = s1_dx + x_dx - z1_dx;
        s1_dy = s1_dy + x_dy - z1_dy;            
        s1_dz = s1_dz + x_dz - z1_dz;  
        
        
        rhs1 = Wy + mu_l*(s_l+real(ifftn(kernel.*fftn(y_l))));
        rhs2 = Wy + mu_e*(s_e+real(ifftn(kernel.*fftn(y_e))))/beta;
        
        z_e = (rhs2)./ a3;
        z_l = (rhs1-weight.*z_e)./ a1;
        %z_l = (a3.*rhs1 - a2.*rhs2)./detA;
        %z_e = (a1.*rhs2 - a2.*rhs1)./detA;
        
        s_e = s_e + real(ifftn(kernel.*fftn(y_e))) - z_e;
        s_l = s_l + real(ifftn(kernel.*fftn(y_l))) - z_l;
        
        
        
        
    end
end
out.time = toc;toc

out.x = mask.*x;
out.phi = real(ifftn(kernel.*fftn(maskc.*x)));
out.iter = t;



end
