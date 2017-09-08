% This is a sample script to show how to use the functions in this toolbox, 
% and how to set the principal variables.
% This example uses an analytic brain phantom.
%
% Based on the code by Bilgic Berkin at http://martinos.org/~berkin/software.html
% Last modified by Carlos Milovic in 2017.03.30
%


%% load data

load chi_phantom
load mask_phantom
load spatial_res

N = size(chi);

imagesc3d2(chi - (mask_use==0), N/2, 1, [90,90,90], [-0.12,0.12], 0, 'True Susceptibility') 
 
center = N/2 + 1;


% Simulate magnitude data

mag = chi-min(chi(:));
mag = mag/max(mag(:));


% Add simulated lesions - Constant spheres

center1 = N/2;
[chiS, Bnuc] = chi_intsphere(spatial_res.*N,N,center1,6, -0.3, 0);
chi = chi+chiS;
mag = mag.*(0.0+1+chiS/0.3);

center2 = [59 86 75];
[chiS, Bnuc] = chi_intsphere(spatial_res.*N,N,center2,6, +0.6, 0);
chi = chi+chiS;
mag = mag.*(0.0+1-chiS/0.6);

center3 = [163 94 24];
[chiS, Bnuc] = chi_intsphere(spatial_res.*N,N,center3,6, +1.2, 0);
chi = chi+chiS;
mag = mag.*(0.0+1-chiS/1.2);

center4 = [105 171 57];
[chiS, Bnuc] = chi_intsphere(spatial_res.*N,N,center4,6, -0.5, 0);
chi = chi+chiS;
mag = mag.*(0.0+1+chiS/0.5);



%% Create dipole kernel and susceptibility to field model

kernel = dipole_kernel( N, spatial_res, 0 ); % 0 for the continuous kernel by Salomir and Marques and Bowtell.
                                             % 1 for a discrete kernel
                                             % formulation
                                             % 2 for an integrated Green
                                             % function

chi = chi-mean(chi(:));
phase_true = ifftn(kernel .* fftn(chi));


%% Add noise

matrix_size = size(chi);
SNR = 345; % peak SNR value
noise = 1/SNR;
signal = mag.*exp(1i*phase_true)+noise*(randn(matrix_size)+1i*randn(matrix_size));
phase_use = angle(signal);
magn_use = abs(signal);


rmse_noise = 100 * norm(mask_use(:) .* (phase_use(:) - phase_true(:))) / norm(mask_use(:).*phase_true(:));


% Add phase 2pi jumps into the phase data to simulare unwrapping errors

pt = [152 102 center(3)];
phase_use(pt(1),pt(2),pt(3)) = phase_use(pt(1),pt(2),pt(3)) - 2*pi;
pt = [110 152 center(3)];
phase_use(pt(1),pt(2),pt(3)) = phase_use(pt(1),pt(2),pt(3)) + 2*pi;
pt = [113 118 center(3)];
phase_use(pt(1),pt(2),pt(3)) = phase_use(pt(1),pt(2),pt(3)) + 2*pi;
pt = [78 106 center(3)];
phase_use(pt(1),pt(2),pt(3)) = phase_use(pt(1),pt(2),pt(3)) - 2*pi;
pt = [133 183 center(3)];
phase_use(pt(1),pt(2),pt(3)) = phase_use(pt(1),pt(2),pt(3)) + 2*pi;

imagesc3d2(phase_use_old, N/2, 2, [90,90,90], [-0.4 0.4], 0, ['Noise RMSE: ', num2str(rmse_noise)]) %-0.04,0.04
imagesc3d2(magn_use, N/2, 3, [90,90,90], [0,1], 0, ['Noisy Magnitude']) 
 



%% TKD recon 

kthre = 0.08;       % truncation threshold

chi_tkd = tkd( phase_use, mask_use, kernel, kthre, N );

rmse_tkd = 100 * norm(real(chi_tkd(:)).*mask_use(:) - chi(:)) / norm(chi(:));
metrics_tkd = compute_metrics(real(chi_tkd.*mask_use),chi);

imagesc3d2(chi_tkd .* mask_use - (mask_use==0), N/2, 3, [90,90,90], [-0.12,0.12], 0, ['TKD RMSE: ', num2str(rmse_tkd)])


%% Closed-form L2 recon

beta = 3e-3;    % regularization parameter

chi_L2 = chiL2( phase_use, mask_use, kernel, beta, N );

rmse_L2 = 100 * norm(real(chi_L2(:)).*mask_use(:) - chi(:)) / norm(chi(:));
metrics_l2 = compute_metrics(real(chi_L2.*mask_use),chi);
imagesc3d2(chi_L2 .* mask_use - (mask_use==0), N/2, 4, [90,90,90], [-0.12,0.12], 0, ['L2 RMSE: ', num2str(rmse_L2)])



%% Linear TV and TGV ADMM recon
% Optional parameters are commented

%num_iter = 50;
%tol_update = 1;

% Common parameters
params = [];

%params.maxOuterIter = num_iter;
%params.tol_update = tol_update;
params.K = kernel;
params.input = phase_use;

%params.mu2 = 1.0;
params.mu1 = 1e-2;                  % gradient consistency
params.alpha1 = 2e-4;               % gradient L1 penalty

% TV
%params.weight = ones(N); 
out = wTV(params); 
rmse_tv = 100 * norm(out.x(:).*mask_use(:) - chi(:)) / norm(chi(:));
metrics_tv = compute_metrics(real(out.x.*mask_use),chi);

params.weight = magn_use;
outw = wTV(params); 
rmse_tvw = 100 * norm(outw.x(:).*mask_use(:) - chi(:)) / norm(chi(:));
metrics_tvw = compute_metrics(real(outw.x.*mask_use),chi);

imagesc3d2(real(out.x) .* mask_use - (mask_use==0), N/2, 5, [90,90,90], [-0.12,0.12], 0, ['TV RMSE: ', num2str(rmse_tv), '  iter : ', num2str(out.iter)])

imagesc3d2(real(outw.x) .* mask_use - (mask_use==0), N/2, 6, [90,90,90], [-0.12,0.12], 0, ['Weighted TV RMSE: ', num2str(rmse_tvw), '  iter : ', num2str(outw.iter)])


% TGV
%params.mu0 = 2*params.mu1;            % second order gradient consistency
%params.alpha0 = 2 * params.alpha1;  % second order gradient L1 penalty

params.weight = ones(N);
out2 = wTGV(params); 
rmse_tgv = 100 * norm(out2.x(:).*mask_use(:) - chi(:)) / norm(chi(:));
metrics_tgv = compute_metrics(real(out2.x.*mask_use),chi);

params.weight = magn_use;
out2w = wTGV(params); 
rmse_tgvw = 100 * norm(out2w.x(:).*mask_use(:) - chi(:)) / norm(chi(:));
metrics_tgvw = compute_metrics(real(out2w.x.*mask_use),chi);

imagesc3d2(real(out2.x) .* mask_use - (mask_use==0), N/2, 7, [90,90,90], [-0.12,0.12], 0, ['TGV RMSE: ', num2str(rmse_tgv), '  iter : ', num2str(out2.iter)])

imagesc3d2(real(out2w.x) .* mask_use - (mask_use==0), N/2, 8, [90,90,90], [-0.12,0.12], 0, ['Weighted TGV RMSE: ', num2str(rmse_tgvw), '  iter : ', num2str(out2w.iter)])




%% Non linear TV and TGV ADMM recon
% Optional parameters are commented
% num_iter = 50;
% tol_update = 1;

% Common parameters
params = [];
% 
% params.maxOuterIter = num_iter;
% params.tol_update = tol_update;
params.K = kernel;
params.input = phase_use;
params.weight = magn_use; % Always required for nonlinear algorithms
 
% params.mu2 = 1.0;
params.mu1 = 1e-2;                  % gradient consistency
params.alpha1 = 2e-4;               % gradient L1 penalty

% nlTV

outnl = nlTV(params); 
rmse_tvnl = 100 * norm(outnl.x(:).*mask_use(:) - chi(:)) / norm(chi(:));
metrics_nltv = compute_metrics(real(outnl.x.*mask_use),chi);

imagesc3d2(real(outnl.x) .* mask_use - (mask_use==0), N/2, 9, [90,90,90], [-0.12,0.12], 0, ['nlTV RMSE: ', num2str(rmse_tvnl), '  iter : ', num2str(outnl.iter)])


% nlTGV
%params.mu0 = 2*params.mu1;            % second order gradient consistency
%params.alpha0 = 2 * params.alpha1;  % second order gradient L1 penalty

out2nl = nlTGV(params); 
rmse_tgvnl = 100 * norm(out2nl.x(:).*mask_use(:) - chi(:)) / norm(chi(:));
metrics_nltgv = compute_metrics(real(out2nl.x.*mask_use),chi);

imagesc3d2(real(out2nl.x) .* mask_use - (mask_use==0), N/2, 10, [90,90,90], [-0.12,0.12], 0, ['nlTV RMSE: ', num2str(rmse_tgvnl), '  iter : ', num2str(out2nl.iter)])


%% Spatially weighted regularization term (example with linear function) - Optional

% Common parameters
params = [];

params.K = kernel;
params.input = phase_use;

params.mu1 = 1e-2;                  % gradient consistency
params.alpha1 = 2e-4;               % gradient L1 penalty

Gm = gradient_calc(magn_use,0); % 0 for vectorized gradient. Use options 1 or 2 to use the L1 or L2 of the gradient.

Gm = max(Gm,noise); % For continous weighting. Soft-threshold to avoid zero divisions and control noise.
params.regweight = mean(Gm(:))./Gm;

outrw = wTV(params); 
rmse_tvrw = 100 * norm(outrw.x(:).*mask_use(:) - chi(:)) / norm(chi(:));
metrics_tvrw = compute_metrics(real(outrw.x.*mask_use),chi);


bGm = threshold_gradient( Gm, mask_use, noise, 0.3 ); % For binary weighting. 
params.regweight = bGm;

outbrw = wTV(params); 
rmse_tvbrw = 100 * norm(outbrw.x(:).*mask_use(:) - chi(:)) / norm(chi(:));
metrics_tvbrw = compute_metrics(real(outbrw.x.*mask_use),chi);

%% FANSI main function call example

options = [];
mu1 = 1e-2;                  % gradient consistency
alpha1 = 2e-4;               % gradient L1 penalty
outf = FANSI( phase_use, magn_use, spatial_res, alpha1, mu1, noise );
rmse_f1 = 100 * norm(outf.x(:).*mask_use(:) - chi(:)) / norm(chi(:));
 [ data_cost1, reg_cost1 ] = compute_costs( outf.x, phase_use, kernel );

options.tgv = true;
options.nonlinear = true;
outf2 = FANSI( phase_use, magn_use, spatial_res, alpha1, mu1, noise, options );
rmse_f2 = 100 * norm(outf2.x(:).*mask_use(:) - chi(:)) / norm(chi(:));
 [ data_cost2, reg_cost2 ] = compute_costs( outf2.x, phase_use, kernel );
