% QSM with a Total Variation regularization with spatially variable fidelity and regularization weights.
% This uses ADMM to solve the functional.
%
% Based on the code by Bilgic Berkin at http://martinos.org/~berkin/software.html
% Last modified by Carlos Milovic in 2017.03.30
%
function out = pdf2wTV(params)
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

    if isfield(params,'mul')
         mul = params.mul;
    else
        mul = 1.0;
    end
    
    if isfield(params,'mue')
         mue = params.mue;
    else
        mue = 1.0;
    end
    
    if isfield(params,'beta')
         beta = params.beta;
    else
        beta = 1.0;
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
    
    if isfield(params,'mask')
        mask = params.mask;
    else
        mask = weight > 0;
    end
    
    if isfield(params,'regweight')
        regweight = params.regweight;
        if length(size(regweight)) == 3
            regweight = repmat(regweight,[1,1,1,3]);
        end
    else
        regweight = ones([N 3]);
    end
    
    
    Wy = weight.*params.input;
    
    
    
    
z_dx = zeros(N, 'single');
z_dy = zeros(N, 'single');
z_dz = zeros(N, 'single');

s_dx = zeros(N, 'single');
s_dy = zeros(N, 'single');
s_dz = zeros(N, 'single');

x = zeros(N, 'single');
xe = zeros(N, 'single');
xl = zeros(N, 'single');

ze = zeros(N, 'single');
zl = zeros(N, 'single');
se = zeros(N, 'single');
sl = zeros(N, 'single');


%     if isfield(params,'precond')
%         precond = params.precond;
%     else
%         precond = true;
%     end
%     
%     if precond
%         z2 = Wy;
%     else
%         z2 = zeros(N,'single');
%     end
% s2 = zeros(N,'single');

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


bx = mu1*EE2;
ax = mue*K2+bx;
dx = mul*K2+bx;
detx = ax.*dx-bx.*bx+eps;

az = weight+mul;
bz = weight;
dz = (1+beta)*weight+mue;
detz = az.*dz-bz.*bz+eps;



% 
% [ky,kx,kz] = meshgrid((-N(2)/2):(N(2)/2-1), (-N(1)/2):(N(1)/2-1), (-N(3)/2):(N(3)/2-1));
% G = (kx).^2+(ky).^2+(kz).^2;
% G = fftshift( exp( -G/(5^2) ) );
% G = G/sum(G(:));
% ze = real(ifftn(fftn(G).*fftn(params.input.*mask)));

%tic
for t = 1:num_iter
    % update x : susceptibility estimate
    tx = E1t .* fftn(z_dx - s_dx);
    ty = E2t .* fftn(z_dy - s_dy);
    tz = E3t .* fftn(z_dz - s_dz);
    T = mu1 * (tx + ty + tz);
    
    x_prev = x;
    De = mue*conj(kernel) .* fftn(ze-se);
    Dl = mul*conj(kernel) .* fftn(zl-sl);
    xe = real(ifftn( (dx.*(De+T)-bx.*(Dl+T))./detx   ));
    xl = real(ifftn( (-bx.*(De+T)+ax.*(Dl+T))./detx   ));
    x = xe+xl;
    xe = (1-mask).*x;
    xl = mask.*x;

    x_update = 100 * norm(x(:)-x_prev(:)) / norm(x(:));
    disp(['Iter: ', num2str(t), '   Update: ', num2str(x_update)])
    
    if x_update < tol_update
        break
    end
    
        ll = alpha1/mu1;
    if t < num_iter
        % update z : gradient varible
        Fx = fftn(x);
        x_dx = real(ifftn(E1 .* Fx));
        x_dy = real(ifftn(E2 .* Fx));
        x_dz = real(ifftn(E3 .* Fx));
        
        z_dx = max(abs(x_dx + s_dx) - regweight(:,:,:,1)*ll, 0) .* sign(x_dx + s_dx);
        z_dy = max(abs(x_dy + s_dy) - regweight(:,:,:,2)*ll, 0) .* sign(x_dy + s_dy);
        z_dz = max(abs(x_dz + s_dz) - regweight(:,:,:,3)*ll, 0) .* sign(x_dz + s_dz);
    
        % update s : Lagrange multiplier
        s_dx = s_dx + x_dx - z_dx;
        s_dy = s_dy + x_dy - z_dy;            
        s_dz = s_dz + x_dz - z_dz;  
        
        
        Pl = real(ifftn( kernel.*fftn(xl)));
        Pe = real(ifftn( kernel.*fftn(xe)));
        ze = ( dz.*(Wy+mul*(Pl+sl))-bz.*((1+beta)*Wy+mue*(Pe+se)) )./detz;
        zl = ( -bz.*(Wy+mul*(Pl+sl))+az.*((1+beta)*Wy+mue*(Pe+se)) )./detz;
        
        se = se + Pe - ze;
        sl = sl + Pl - zl;
    end
end
out.time = toc;toc

out.x = xl;
out.phi = Pe;
out.iter = t;



end
