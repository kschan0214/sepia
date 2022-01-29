function out = WH_nlTV_4sepia(params)
% Nonlinear Weak Harmonics - QSM and Total Variation regularization 
% with spatially variable fidelity and regularization weights.
% This uses ADMM to solve the functional.
% This function is used to remove background field remnants from *local* field maps and
% calculate the susceptibility of tissues simultaneously.
%
% Parameters: params - structure with 
% Required fields:
% params.input: local field map
% params.K: dipole kernel in the frequency space
% params.alpha1: gradient penalty (L1-norm) or regularization weight
% Optional fields:
% params.beta: harmonic constrain weight (default value = 150)
% params.muh: harmonic consistency weight (recommended value = beta/50)
% params.mask: ROI to calculate susceptibility values (if not provided, will be calculated from 'weight')
% params.mu1: gradient consistency weight (ADMM weight, recommended = 100*alpha1)
% params.mu2: fidelity consistency weight (ADMM weight, recommended value = 1.0)
% params.maxOuterIter: maximum number of iterations (recommended for testing = 150, for correct 
%                      convergence of the harmonic field hundreds of iterations are needed)
% params.tol_update: convergence limit, update ratio of the solution (recommended = 0.1)
% params.weight: data fidelity spatially variable weight (recommended = magnitude_data). 
% params.regweight: regularization spatially variable weight.
% params.precond: preconditionate solution by smart initialization
%
% Output: out - structure with the following fields:
% out.x: calculated susceptibility map
% out.phi: harmonic phase in [same range as input]
% out.iter: number of iterations needed
% out.time: total elapsed time (including pre-calculations)
%
% Modified by Carlos Milovic in 2017.06.09
% Last modified by Carlos Milovic in 2020.07.11
%

tic

    % Required parameters
alpha = params.alpha1;
kernel = params.K;
phase = params.input;

N = size(params.input);
    
    % Optional parameters
    if isfield(params,'mu1')
         mu = params.mu1;
    else
        mu = 100*alpha;
    end
    if isfield(params,'mu2')
         mu2 = params.mu2;
    else
        mu2 = 1.0;
    end
    
    if isfield(params,'muh')
         muh = params.muh;
    else
        muh = 5.0;
    end
    if isfield(params,'beta')
         beta = params.beta;
    else
        beta = 150.0;
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
        mask =weight > 0;
    end
    

    if isfield(params,'maxOuterIter')
        num_iter = params.maxOuterIter;
    else
        num_iter = 150;
    end
    
    if isfield(params,'tol_update')
       tol_update  = params.tol_update;
    else
       tol_update = 0.1;
    end

    
    if isfield(params,'regweight')
        regweight = params.regweight;
        if length(size(regweight)) == 3
            regweight = repmat(regweight,[1,1,1,3]);
        end
    else
        regweight = ones([N 3]);
    end
    
    if isfield(params,'delta_tol')
        delta_tol = params.delta_tol;
    else
        delta_tol = 1e-6;
        
    end
    
    if isfield(params,'precond')
        precond = params.precond;
    else
        precond = true;
    end
    
    
    % Variable initialization
z_dx = zeros(N, 'single');
z_dy = zeros(N, 'single');
z_dz = zeros(N, 'single');

s_dx = zeros(N, 'single');
s_dy = zeros(N, 'single');
s_dz = zeros(N, 'single');

x = zeros(N, 'single');

phi_h = zeros(N, 'single');

z_h = zeros(N, 'single');
s_h = zeros(N, 'single');
% KC 20211114, Wy is unknown to me yet
%     if precond
%         z2 = Wy;
%     else
        z2 = zeros(N,'single');
%     end
s2 = zeros(N,'single');

kernel = params.K;


% Define the operators
[k1, k2, k3] = ndgrid(0:N(1)-1,0:N(2)-1,0:N(3)-1);

E1 = 1 - exp(2i .* pi .* k1 / N(1));
E2 = 1 - exp(2i .* pi .* k2 / N(2));
E3 = 1 - exp(2i .* pi .* k3 / N(3));

E1t = conj(E1);
E2t = conj(E2);
E3t = conj(E3);

EE2 = E1t .* E1 + E2t .* E2 + E3t .* E3;
Lap = EE2;%E1+E1t+E2+E2t+E3+E3t;
K2 = abs(kernel).^2;

%tic
ll = alpha/mu;
for t = 1:num_iter
    
   
    % update x : susceptibility estimate
    tx = E1t .* fftn(z_dx - s_dx);
    ty = E2t .* fftn(z_dy - s_dy);
    tz = E3t .* fftn(z_dz - s_dz);
    
    x_prev = x;
    Dt_kspace = conj(kernel) .* fftn(z2-s2-(phi_h));
    x = mask.*real(ifftn( (mu * (tx + ty + tz) + Dt_kspace) ./ (eps + mu2*K2 + mu * EE2) ));

    x_update = 100 * norm(x(:)-x_prev(:)) / norm(x(:));
    disp(['Iter: ', num2str(t), '   Update: ', num2str(x_update)])
    
    if x_update < tol_update
        break
    end
    
    
    if t < num_iter
        % update z : gradient variable
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
        
        
        
        rhs_z2 = mu2*real(ifftn(kernel.*Fx)+s2 +phi_h  );
        z2 =  rhs_z2 ./ mu2 ;

        % Newton-Raphson method
        delta = inf;
        inn = 0;
        while (delta > delta_tol && inn < 50)
            inn = inn + 1;
            norm_old = norm(z2(:));
            
            update = ( weight .* sin(z2 - phase) + mu2*z2 - rhs_z2 ) ./ ( weight .* cos(z2 - phase) + mu2 );            
        
            z2 = z2 - update;     
            delta = norm(update(:)) / norm_old;
        end        
        %disp(delta)
    disp(['Inner Solver - Iter: ', num2str(inn), '   Update: ', num2str(delta)])
        
               
        
        Fphi_h = (muh * conj(Lap).*fftn(z_h-s_h) + mu2*fftn(z2-s2) - mu2*kernel.*Fx) ./ (eps + mu2 + muh * Lap.*conj(Lap)) ;
        phi_h = real(ifftn(Fphi_h));
        
        z_h = muh*(real(ifftn(Lap.*Fphi_h))+s_h)./(muh+beta*mask);
        
        s2 = s2 + real(ifftn(kernel.*Fx)) - z2+phi_h;
        s_h = s_h + real(ifftn(Lap.*Fphi_h)) - z_h;
        
    end
    
    
end
% Extract output values
out.time = toc;toc

out.x = x;
out.phi = phi_h;
out.iter = t;



end
% % Single-Step QSM with a Total Variation regularization with spatially variable fidelity and regularization weights.
% % This uses ADMM to solve the functional.
% %
% % Last modified by Carlos Milovic in 2017.06.09
% % bug fix: line 165 (KC, 2 March 2020)
% %
% function out = WH_nlTV_4sepia(params)
% %
% % input: params - structure with the following required fields:
% % params.input - local field map
% % params.K - dipole kernel in the frequency space
% % params.alpha1 - gradient L1 penalty, regularization weight
% % params.mu1 - gradient consistency weight
% %                 and the following optional fields:
% % params.mu2 - fidelity consistency weight (recommended value = 1.0)
% % params.beta - harmonic constrain weight
% % params.muh - harmonic consistency weight (recommended value = beta)
% % params.mask - ROI to calculate susceptibility values (if not provided, will be calculated from 'weight')
% % params.maxOuterIter - maximum number of ADMM iterations (recommended = 50)
% % params.tol_update - convergence limit, change rate in the solution (recommended = 1.0)
% % params.weight - data fidelity spatially variable weight (recommended = magnitude_data). Not used if not specified
% % params.regweight - regularization spatially variable weight. Not used if not specified
% % params.N - array size
% % params.precond - preconditionate solution (for stability)
% %
% % output: out - structure with the following fields:
% % out.x - calculated susceptibility map
% % out.iter - number of iterations needed
% % out.time - total elapsed time (including pre-calculations)
% % out.phi - harmonic phase in radians
% %
% 
% 
% tic
% 
% mu = params.mu1;
% alpha = params.alpha1;
% 
%     if isfield(params,'mu2')
%          mu2 = params.mu2;
%     else
%         mu2 = 1.0;
%     end
%     
%     if isfield(params,'muh')
%          muh = params.muh;
%     else
%         muh = 1.0;
%     end
%     if isfield(params,'beta')
%          beta = params.beta;
%     else
%         beta = 1.0;
%     end
%     
%     if isfield(params,'mask')
%          mask = params.mask;
%     else
%         mask = params.weight > 0;
%     end
%     
%     if isfield(params,'N')
%          N = params.N;
%     else
%         N = size(params.input);
%     end
% 
%     if isfield(params,'maxOuterIter')
%         num_iter = params.maxOuterIter;
%     else
%         num_iter = 50;
%     end
%     
%     if isfield(params,'tol_update')
%        tol_update  = params.tol_update;
%     else
%        tol_update = 1;
%     end
% 
%     if isfield(params,'weight')
%         weight = params.weight;
%     else
%         weight = ones(N);
%     end
%     weight = weight.*weight;
%     
%     if isfield(params,'regweight')
%         regweight = params.regweight;
%         if length(size(regweight)) == 3
%             regweight = repmat(regweight,[1,1,1,3]);
%         end
%     else
%         regweight = ones([N 3]);
%     end
%     
%     if ~isfield(params,'delta_tol')
%         delta_tol = 1e-6;
%     else
%         delta_tol = params.delta_tol;
%     end
%     
%     %Wy = (weight.*params.input./(weight+mu2));
%     
% z_dx = zeros(N, 'single');
% z_dy = zeros(N, 'single');
% z_dz = zeros(N, 'single');
% 
% s_dx = zeros(N, 'single');
% s_dy = zeros(N, 'single');
% s_dz = zeros(N, 'single');
% 
% x = zeros(N, 'single');
% 
% phi_h = zeros(N, 'single');
% 
% z_h = zeros(N, 'single');
% s_h = zeros(N, 'single');
% % 
% %     if isfield(params,'precond')
% %         precond = params.precond;
% %     else
% %         precond = true;
% %     end
% %     
% %     if precond
% %         z2 = Wy;
% %     else
% %         z2 = zeros(N,'single');
% %     end
% %         z2 = zeros(N,'single');
% s2 = zeros(N,'single');
% 
% kernel = params.K;
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
% Lap = E1+E1t+E2+E2t+E3+E3t;
% K2 = abs(kernel).^2;
% % 
% % phi_e = real(ifftn( kernel.*fftn(1-mask) ));
% % lambda = 0;
% % phi0 = 0;
% % 
% %     a1 = sum( mask(:).*phi_e(:) );
% %     a2 = sum( mask(:).*phi_e(:).*phi_e(:) );
% %     a3 = sum( mask(:) );
% %     f1 = sum( mask(:).*phi_e(:).*params.input(:) );
% %     f2 = sum( mask(:).*params.input(:) );
% %     det = a1*a1-a2*a3;
% %     
% %     phi0 = (a1*f1 - a2*f2)/det;
% %     lambda = (-a3*f1 + a1*f2)/det;
% %     PHI = lambda*phi_e+phi0;
% 
%     if isfield(params,'background') %if ~isfield(params,'background')
%         PHI = params.background;
%     else
%         PHI = zeros(N);
%     end
% %z2 = PHI;
% phase = (params.input-PHI);
%     %Wy = (weight.*(params.input-PHI)./(weight+mu2));
% %z2 = weight.*(unwrapLaplacian(params.input.*mask,N)-PHI);
% z2 = weight.*params.input;
% %z2 = zeros(N,'single');
% %tic
% for t = 1:num_iter
%     
%    
%     % update x : susceptibility estimate
%     tx = E1t .* fftn(z_dx - s_dx);
%     ty = E2t .* fftn(z_dy - s_dy);
%     tz = E3t .* fftn(z_dz - s_dz);
%     
%     x_prev = x;
%     Dt_kspace = conj(kernel) .* fftn(z2-s2-(phi_h));
%     x = mask.*real(ifftn( (mu * (tx + ty + tz) + Dt_kspace) ./ (eps + mu2*K2 + mu * EE2) ));
% 
%     x_update = 100 * norm(x(:)-x_prev(:)) / norm(x(:));
%     disp(['Iter: ', num2str(t), '   Update: ', num2str(x_update)])
%     
%     if x_update < tol_update
%         break
%     end
%     
%     
%         ll = alpha/mu;
%     if t < num_iter
%         % update z : gradient varible
%         Fx = fftn(x);
%         x_dx = real(ifftn(E1 .* Fx));
%         x_dy = real(ifftn(E2 .* Fx));
%         x_dz = real(ifftn(E3 .* Fx));
%         
%         z_dx = max(abs(x_dx + s_dx) - regweight(:,:,:,1)*ll, 0) .* sign(x_dx + s_dx);
%         z_dy = max(abs(x_dy + s_dy) - regweight(:,:,:,2)*ll, 0) .* sign(x_dy + s_dy);
%         z_dz = max(abs(x_dz + s_dz) - regweight(:,:,:,3)*ll, 0) .* sign(x_dz + s_dz);
%     
%         % update s : Lagrange multiplier
%         s_dx = s_dx + x_dx - z_dx;
%         s_dy = s_dy + x_dy - z_dy;            
%         s_dz = s_dz + x_dz - z_dz;  
%         
%         
%         
%         rhs_z2 = mu2*real(ifftn(kernel.*Fx)+s2 +phi_h  );
%         z2 =  rhs_z2 ./ mu2 ;
% 
%         % Newton-Raphson method
%         delta = inf;
%         inn = 0;
%         while (delta > delta_tol && inn < 50)
%             inn = inn + 1;
%             norm_old = norm(z2(:));
%             
%             update = ( weight .* sin(z2 - phase) + mu2*z2 - rhs_z2 ) ./ ( weight .* cos(z2 - phase) + mu2 );            
%         
%             z2 = z2 - update;     
%             delta = norm(update(:)) / norm_old;
%         end        
%         %disp(delta)
%     disp(['Inner Solver - Iter: ', num2str(inn), '   Update: ', num2str(delta)])
%         
%         %z2 = Wy + mu2*real(ifftn(kernel.*Fx)+s2+phi_h )./(weight + mu2);
%                
%         
%         Fphi_h = (muh * conj(Lap).*fftn(z_h-s_h) + mu2*fftn(z2-s2) - mu2*kernel.*Fx) ./ (eps + mu2 + muh * Lap.*conj(Lap)) ;
%         phi_h = real(ifftn(Fphi_h));
%         %phi_h = phi_h-mean( phi_h(:) );
%         
%         z_h = muh*(real(ifftn(Lap.*Fphi_h))+s_h)./(muh+beta*mask);
%         
%         s2 = s2 + real(ifftn(kernel.*Fx)) - z2+phi_h;
%         s_h = s_h + real(ifftn(Lap.*Fphi_h)) - z_h;
%         
%         
%         
% %     phi_x = real(ifftn(kernel.*fftn(x)));
% %     
% %     a1 = sum( phi_e(:) );
% %     a2 = sum( phi_e(:).*phi_e(:) );
% %     a3 = numel(x);%sum( mask(:) );
% %     f1 = sum( phi_e(:).*(z2(:)-s2(:)-phi_x(:)-phi_h(:)) );
% %     f2 = sum( (z2(:)-s2(:)-phi_x(:)-phi_h(:)) );
% %     det = a1*a1-a2*a3;
% %     
% %     phi0 = (a1*f1 - a2*f2)/det;
% %     lambda = (-a3*f1 + a1*f2)/det;
% %     PHI = lambda*phi_e+phi0;
%     end
%     
%     
%     
% end
% out.time = toc;toc
% 
% out.x = x;
% out.phi = phi_h;
% out.iter = t;
% 
% 
% 
% end
