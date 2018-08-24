% This is a sample script to show how to use the functions in this toolbox, 
% and how to set the principal variables.
% This example uses an analytic brain phantom.
%
% Based on the code by Bilgic Berkin at http://martinos.org/~berkin/software.html
% Last modified by Carlos Milovic in 2017.03.31
%

set(0,'DefaultFigureWindowStyle','docked')


addpath(genpath(pwd))

%% load data

load chi_phantom
load mask_phantom
load spatial_res_phantom

N = size(chi);
chi_true = mask_use.*(chi-( sum(mask_use(:).*chi(:))/sum(mask_use(:)) ));
imagesc3d2(chi_true, N/2, 1, [90,90,90], [-0.12,0.12], 0, 'True Susceptibility') 
 
center = N/2 + 1;

%chi_true = chi-mean(chi(:));
% Simulate magnitude data
%mag = ones(N);
%mag = chi-min(chi(:));
%mag = mask_use.*mag/max(mag(:));
% 
% chi = chi_true - (1-mask_use)*1.5;
% chi = chi-mean(chi(:));
%Add susceptibility sources outside the ROI - Constant spheres

center1 = [256-32 256-80 49];
[chiS, Bnuc] = chi_intsphere(spatial_res.*N,N,center1,12, -030/10, 0);
chi = chi+chiS;

center2 = [256-233 256-99 49];
[chiS, Bnuc] = chi_intsphere(spatial_res.*N,N,center2,12, +030/10, 0);
chi = chi+chiS;

center3 = [128 168 98-77];
[chiS, Bnuc] = chi_intsphere(spatial_res.*N,N,center3,12, +120/10, 0);
chi = chi+chiS;

center4 = [128 26 98-19];
[chiS, Bnuc] = chi_intsphere(spatial_res.*N,N,center4,12, -120/10, 0);
chi = chi+chiS;

center5 = [228 128 98-62];
[chiS, Bnuc] = chi_intsphere(spatial_res.*N,N,center5,12, 150/10, 0);
chi = chi+chiS;

center6 = [26 128 98-26];
[chiS, Bnuc] = chi_intsphere(spatial_res.*N,N,center6,12, -150/10, 0);
chi = chi+chiS;
imagesc3d2(chi, N/2, 2, [90,90,90], [-0.12,0.12], 0, 'Full Susceptibility') 



xt = zeros( [256+32 256+32 98+30+16] );
xt(17:(256+16),17:(256+16),34:(33+98)) = chi_true;
x = zeros( [256+32 256+32 98+30+16] );
x(17:(256+16),17:(256+16),34:(33+98)) = chi;
N = size(x);
chi_true = xt;
chi = x;
msk = zeros( [256+32 256+32 98+30+16] );
msk(17:(256+16),17:(256+16),34:(33+98)) = mask_use;
mask_use = msk;
mag = chi_true-min(chi_true(:));
mag = mask_use.*mag/max(mag(:));
clear xt x msk;

%% Create dipole kernel and susceptibility to field model

kernel = dipole_kernel( N, spatial_res, 0 ); % 0 for the continuous kernel by Salomir and Marques and Bowtell.
                                             % 1 for a discrete kernel
                                             % formulation
                                             % 2 for an integrated Green
                                             % function

chi = chi-mean(chi(:));
phase_true = ifftn(kernel .* fftn(chi_true));
phase_full = ifftn(kernel .* fftn(chi));

imagesc3d2(phase_true, N/2, 3, [90,90,90], [-0.12,0.12], 0, 'Local Phase') 
imagesc3d2(phase_full, N/2, 4, [90,90,90], [-1.12,1.12], 0, 'Full Phase') 
imagesc3d2((phase_full-phase_true), N/2, 5, [90,90,90], [-1.12,1.12], 0, 'H Phase') 
imagesc3d2((phase_full-phase_true), N/2, 6, [90,90,90], [-0.12,0.12], 0, 'H Phase') 

%% Add noise

matrix_size = size(chi);
SNR = 345; % peak SNR value
noise = 1/SNR;
signal = mag.*exp(1i*(phase_full))+noise*(randn(matrix_size)+1i*randn(matrix_size));
phase_use = angle(signal);
magn_use = abs(signal);


rmse_noise = 100 * norm(mask_use(:) .* (phase_use(:) - phase_full(:))) / norm(mask_use(:).*phase_true(:));


% Add phase 2pi jumps into the phase data to simulate unwrapping errors
% 
% pt = [152 102 center(3)];
% phase_use(pt(1),pt(2),pt(3)) = phase_use(pt(1),pt(2),pt(3)) - 2*pi;
% pt = [110 152 center(3)];
% phase_use(pt(1),pt(2),pt(3)) = phase_use(pt(1),pt(2),pt(3)) + 2*pi;
% pt = [113 118 center(3)];
% phase_use(pt(1),pt(2),pt(3)) = phase_use(pt(1),pt(2),pt(3)) + 2*pi;
% pt = [78 106 center(3)];
% phase_use(pt(1),pt(2),pt(3)) = phase_use(pt(1),pt(2),pt(3)) - 2*pi;
% pt = [133 183 center(3)];
% phase_use(pt(1),pt(2),pt(3)) = phase_use(pt(1),pt(2),pt(3)) + 2*pi;

imagesc3d2(phase_use, N/2, 7, [90,90,90], [-0.12 0.12], 0, ['Noise RMSE: ', num2str(rmse_noise)]) %-0.04,0.04
imagesc3d2(magn_use, N/2, 8, [90,90,90], [0,1], 0, ['Noisy Magnitude']) 
 
imagesc3d2(phase_unwrap, N/2, 8, [90,90,90], [-0.12 0.12], 0, ['Noise RMSE: ', num2str(rmse_noise)]) %-0.04,0.04

imagesc3d2(pu, N/2,9, [90,90,90], [-0.12 0.12], 0, ['Noise RMSE: ', num2str(rmse_noise)]) %-0.04,0.04


imagesc3d2(phase_full, N/2, 7, [90,90,90], [-0.12 0.12], 0, ['Noise RMSE: ', num2str(rmse_noise)]) %-0.04,0.04

imagesc3d2(pt_uw, N/2, 8, [90,90,90], [-0.12 0.12], 0, ['Noise RMSE: ', num2str(rmse_noise)]) %-0.04,0.04

imagesc3d2(pt, N/2,9, [90,90,90], [-0.12 0.12], 0, ['Noise RMSE: ', num2str(rmse_noise)]) %-0.04,0.04

signal2 = mag.*exp(1i*(10*phase_full))+noise*(randn(matrix_size)+1i*randn(matrix_size));
phase_use2 = angle(signal2);

imagesc3d2(phase_use2/10, N/2, 8, [90,90,90], [-0.12 0.12], 0, ['Noise RMSE: ', num2str(rmse_noise)]) %-0.04,0.04
phase_unwrap = unwrapLaplacian(phase_use2, N,spatial_res);
imagesc3d2(phase_unwrap/10, N/2, 9, [90,90,90], [-0.12 0.12], 0, ['Noise RMSE: ', num2str(rmse_noise)]) %-0.04,0.04
phase_unwrap = unwrapLaplacian(phase_use2.*mask_use, N,spatial_res);
imagesc3d2(phase_unwrap/10, N/2, 10, [90,90,90], [-0.12 0.12], 0, ['Noise RMSE: ', num2str(rmse_noise)]) %-0.04,0.04
pu = unwrap(phase_use2, spatial_res);
imagesc3d2(pu/10, N/2,11, [90,90,90], [-0.12 0.12], 0, ['Noise RMSE: ', num2str(rmse_noise)]) %-0.04,0.04


%% TKD recon 

kthre = 0.08;       % truncation threshold

%chi_tkd = tkd( phase_use.*mask_use, mask_use, kernel, kthre, N );
chi_tkd = tkd( phase_use, mask_use, kernel, kthre, N );

rmse_tkd = 100 * norm(real(chi_tkd(:)).*mask_use(:) - chi_true(:)) / norm(chi_true(:));
metrics_tkd = compute_metrics(real(chi_tkd.*mask_use),chi_true);

imagesc3d2(chi_tkd .* mask_use - (mask_use==0), N/2, 9, [90,90,90], [-0.12,0.12], 0, ['TKD RMSE: ', num2str(rmse_tkd)])


%% Closed-form L2 recon

beta = 3e-3;    % regularization parameter

chi_L2 = chiL2( phase_use, mask_use, kernel, beta, N );

rmse_L2 = 100 * norm(real(chi_L2(:)).*mask_use(:) - chi_true(:)) / norm(chi_true(:));
metrics_l2 = compute_metrics(real(chi_L2.*mask_use),chi_true);
imagesc3d2(chi_L2 .* mask_use - (mask_use==0), N/2, 10, [90,90,90], [-0.12,0.12], 0, ['L2 RMSE: ', num2str(rmse_L2)])



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
params.weight = ones(N); 
out = wTV(params); 
rmse_tv = 100 * norm(out.x(:).*mask_use(:) - chi_true(:)) / norm(chi_true(:));
metrics_tv = compute_metrics(real(out.x.*mask_use),chi_true);

imagesc3d2(real(out.x) .* mask_use - (mask_use==0), N/2, 10, [90,90,90], [-0.12,0.12], 0, ['TV RMSE: ', num2str(rmse_tv), '  iter : ', num2str(out.iter)])

params.weight = magn_use;
outw = wTV(params); 
rmse_tvw = 100 * norm(outw.x(:).*mask_use(:) - chi_true(:)) / norm(chi_true(:));
metrics_tvw = compute_metrics(real(outw.x.*mask_use),chi_true);

imagesc3d2(real(outw.x) .* mask_use - (mask_use==0), N/2, 11, [90,90,90], [-0.12,0.12], 0, ['Weighted TV RMSE: ', num2str(rmse_tvw), '  iter : ', num2str(outw.iter)])

save testviernes
% Single step optimization

params.maxOuterIter = 100;
params.tol_update = 0.5;
params.mask = mask_use;
  
%muh = 10.^( ((1:17)-9)/2 );%10.^( ((1:25)-15)/4 +2);
beta = 10.^( ((1:15)-9)/2 );
for sb = 1:15
%for sm = 1:17
%     
 params.beta = beta(sb);
 params.muh = 100*beta(sb);%muh(sm);
%  
% outss = sswTV(params); 
% chiss = outss.x;
% 
% rmse(sb) = compute_rmse(real(chiss.*mask_use),chi_true.*mask_use);
% ssim(sb) = compute_ssim(real(chiss.*mask_use),chi_true.*mask_use);
 
outss = sswTVb(params); 
chiss = outss.x;

rmse2(sb) = compute_rmse(real(chiss.*mask_use),chi_true.*mask_use);
ssim2(sb) = compute_ssim(real(chiss.*mask_use),chi_true.*mask_use);
% 
% outss = sswTVc(params); 
% chiss = outss.x;
% 
% rmse3(sb) = compute_rmse(real(chiss.*mask_use),chi_true.*mask_use);
% ssim3(sb) = compute_ssim(real(chiss.*mask_use),chi_true.*mask_use);
clc
display(sb);
end
% 
% save phantom_ss2
% clc
%end

save phantom_ss2b

[a b] = min(rmse_map4);
[c d] = min(a);
[d b(d) c]

%rmse - 0.8
%ssim - 100

figure;
plot(rmse(1:6),'r');
hold on;
plot(rmse2(1:6),'g');

plot(rmse3(1:6),'b');
hold off;

figure(20)
imshow(rmse_map,[]);



params.beta = beta(8);
params.muh = 100*beta(8);
outss = sswTV(params); 
chiss = outss.x;
rmse_tvss = compute_rmse(real(chiss.*mask_use),chi_true.*mask_use)
imagesc3d2(real(outss.x) .* mask_use - (mask_use==0), N/2, 22, [90,90,90], [-0.12,0.12], 0, ['SS TV RMSE: ', num2str(rmse_tvss), '  iter : ', num2str(outss.iter)])
imagesc3d2(real(outss.phi), N/2, 23, [90,90,90], [-0.12,0.12], 0, ['SS Phi'])

%params.beta = 6;
%params.muh = 1;
params.beta = 30;
params.muh = 100*30;
outss = sswTV(params); 
chiss = outss.x;
rmse_tvssb = compute_rmse(real(chiss.*mask_use),chi_true.*mask_use)
imagesc3d2(real(outss.x) .* mask_use - (mask_use==0), N/2, 12, [90,90,90], [-0.12,0.12], 0, ['SS TV RMSE: ', num2str(rmse_tvssb), '  iter : ', num2str(outss.iter)])
imagesc3d2(real(outss.phi), N/2, 13, [90,90,90], [-0.12,0.12], 0, ['SS Phi'])

for sb = 1:15
%for sm = 1:17
%     
 params.mum = beta(sb);
outss = sswTVm(params); 
chiss = outss.x;

rmsem(sb) = compute_rmse(real(chiss.*mask_use),chi_true.*mask_use);
ssimm(sb) = compute_ssim(real(chiss.*mask_use),chi_true.*mask_use);
% 
% outss = sswTVc(params); 
% chiss = outss.x;
% 
% rmse3(sb) = compute_rmse(real(chiss.*mask_use),chi_true.*mask_use);
% ssim3(sb) = compute_ssim(real(chiss.*mask_use),chi_true.*mask_use);
clc
display(sb);
end


params.gamma = 1;
params.beta = beta(8);
params.muh = 100*beta(8);
params.mu2 = 1;
outss = cwhwTV(params); 
chiss = outss.x;
rmse_tvssc = compute_rmse(real(chiss.*mask_use),chi_true.*mask_use)
imagesc3d2(real(outss.x) .* mask_use - (mask_use==0), N/2, 14, [90,90,90], [-0.12,0.12], 0, ['SS TV RMSE: ', num2str(rmse_tvssc), '  iter : ', num2str(outss.iter)])
imagesc3d2(real(outss.phi), N/2, 15, [90,90,90], [-0.12,0.12], 0, ['SS Phi'])

params.gamma = 1;
params.beta = 100;%beta(8);
params.muh = 1;%beta(8);
params.mu2 = 1;
outss = cwhwTVb(params); 
chiss = outss.x;
rmse_tvssd = compute_rmse(real(chiss.*mask_use),chi_true.*mask_use)
imagesc3d2(real(outss.x) .* mask_use - (mask_use==0), N/2, 16, [90,90,90], [-0.12,0.12], 0, ['SS TV RMSE: ', num2str(rmse_tvssd), '  iter : ', num2str(outss.iter)])
imagesc3d2(real(outss.phi), N/2, 17, [90,90,90], [-0.12,0.12], 0, ['SS Phi'])


params.beta = beta(7);
params.muh = beta(7)*100;
outss = sswTVc(params); 
chiss = outss.x;
rmse_tvss = compute_rmse(real(chiss.*mask_use),chi_true.*mask_use);
imagesc3d2(real(outss.x) .* mask_use - (mask_use==0), N/2, 14, [90,90,90], [-0.12,0.12], 0, ['SS TV RMSE: ', num2str(rmse_tvss), '  iter : ', num2str(outss.iter)])

imagesc3d2(real(outss.phi), N/2, 15, [90,90,90], [-0.12,0.12], 0, ['SS Phi'])


params.maxOuterIter = 200;
params.tol_update = 0.1;

params.input = phase_use;
params.weight = magn_use.*mask_use;

params.beta = (1)*1;
params.mu_e = (1)*1;
params.mu_l = 1;
params.nu_e = 1/300;
params.nu_l = 1/300;
outp = wPDFssTV( params );
chiss = outp.x;
rmse_tvp = compute_rmse(real(chiss.*mask_use),chi_true.*mask_use);
imagesc3d2(real(outp.x) .* mask_use - (mask_use==0), N/2, 31, [90,90,90], [-0.12,0.12], 0, ['P TV RMSE: ', num2str(rmse_tvp), '  iter : ', num2str(outp.iter)])

imagesc3d2(real(outp.phi), N/2, 32, [90,90,90], [-0.12,0.12], 0, ['SS Phi'])

outp = wPDFssTVb( params );
chiss = outp.x;
rmse_tvp = compute_rmse(real(chiss.*mask_use),chi_true.*mask_use);
imagesc3d2(real(outp.x) .* mask_use - (mask_use==0), N/2, 33, [90,90,90], [-0.12,0.12], 0, ['P TV RMSE: ', num2str(rmse_tvp), '  iter : ', num2str(outp.iter)])

imagesc3d2(real(outp.phi), N/2, 34, [90,90,90], [-0.12,0.12], 0, ['SS Phi'])
outp = wPDFssTVc( params );
chiss = outp.x;
rmse_tvp = compute_rmse(real(chiss.*mask_use),chi_true.*mask_use);
imagesc3d2(real(outp.x) .* mask_use - (mask_use==0), N/2, 35, [90,90,90], [-0.12,0.12], 0, ['P TV RMSE: ', num2str(rmse_tvp), '  iter : ', num2str(outp.iter)])

imagesc3d2(real(outp.phi), N/2, 36, [90,90,90], [-0.12,0.12], 0, ['SS Phi'])
%%
params.maxOuterIter = 500;
params.tol_update = 0.1;
outss = sswTV(params); 
chiss = outss.x;
rmse_tvss = compute_rmse(real(chiss.*mask_use),chi_true.*mask_use);
imagesc3d2(real(outss.x) .* mask_use - (mask_use==0), N/2, 16, [90,90,90], [-0.12,0.12], 0, ['SS TV RMSE: ', num2str(rmse_tvss), '  iter : ', num2str(outss.iter)])

imagesc3d2(real(outss.phi), N/2, 17, [90,90,90], [-0.12,0.12], 0, ['SS Phi'])

params.beta = 1;
params.muh = 1;
outss = sswTV(params); 
chiss = outss.x;
rmse_tvss = compute_rmse(real(chiss.*mask_use),chi_true.*mask_use);
imagesc3d2(real(outss.x) .* mask_use - (mask_use==0), N/2, 18, [90,90,90], [-0.12,0.12], 0, ['SS TV RMSE: ', num2str(rmse_tvss), '  iter : ', num2str(outss.iter)])

imagesc3d2(real(outss.phi), N/2, 19, [90,90,90], [-0.12,0.12], 0, ['SS Phi'])

params.maxOuterIter = 200;
params.tol_update = 0.1;
outss = sswTV(params); 
chiss = outss.x;
rmse_tvss = compute_rmse(real(chiss.*mask_use),chi_true.*mask_use);
imagesc3d2(real(outss.x) .* mask_use - (mask_use==0), N/2, 31, [90,90,90], [-0.12,0.12], 0, ['SS TV RMSE: ', num2str(rmse_tvss), '  iter : ', num2str(outss.iter)])

imagesc3d2(real(outss.phi), N/2, 32, [90,90,90], [-0.12,0.12], 0, ['SS Phi'])


params.beta = beta(1);
params.muh = muh(7);
outss = sswTV(params); 
chiss = outss.x;
rmse_tvss = compute_rmse(real(chiss.*mask_use),chi_true.*mask_use);
imagesc3d2(real(outss.x) .* mask_use - (mask_use==0), N/2, 33, [90,90,90], [-0.12,0.12], 0, ['SS TV RMSE: ', num2str(rmse_tvss), '  iter : ', num2str(outss.iter)])

imagesc3d2(real(outss.phi), N/2, 34, [90,90,90], [-0.12,0.12], 0, ['SS Phi'])

params.maxOuterIter = 200;
params.tol_update = 0.1;

params.beta = beta(4);
params.muh = 2*muh(8);
outss = sswTV(params); 
chiss = outss.x;
rmse_tvss = compute_rmse(real(chiss.*mask_use),chi_true.*mask_use);
imagesc3d2(real(outss.x) .* mask_use - (mask_use==0), N/2, 35, [90,90,90], [-0.12,0.12], 0, ['SS TV RMSE: ', num2str(rmse_tvss), '  iter : ', num2str(outss.iter)])

imagesc3d2(real(outss.phi), N/2, 36, [90,90,90], [-0.12,0.12], 0, ['SS Phi'])

params.mask = mask_use;
outssb = efwTV(params); 
chissb = outssb.x;
rmse_tvssb = compute_rmse(real(chissb.*mask_use),chi_true.*mask_use);
imagesc3d2(real(outssb.x) .* mask_use - (mask_use==0), N/2, 37, [90,90,90], [-1.12,1.12], 0, ['SS TVb RMSE: ', num2str(rmse_tvssb), '  iter : ', num2str(outssb.iter)])

imagesc3d2(real(outssb.phi), N/2, 38, [90,90,90], [-1.12,1.12], 0, ['SS Phib'])

params.mask = mask_use;
outssb = sswTVb(params); 
chissb = outssb.x;
rmse_tvssb = compute_rmse(real(chissb.*mask_use),chi_true.*mask_use);
imagesc3d2(real(outssb.x) .* mask_use - (mask_use==0), N/2, 37, [90,90,90], [-1.12,1.12], 0, ['SS TVb RMSE: ', num2str(rmse_tvssb), '  iter : ', num2str(outssb.iter)])

imagesc3d2(real(outssb.phi), N/2, 38, [90,90,90], [-1.12,1.12], 0, ['SS Phib'])


Lpf = real(ifftn( Lap.*fftn(phase_full) ));
Lpx = real(ifftn( Lap.*fftn(phase_true) ));
Lpu = real(ifftn( Lap.*fftn(phase_use) ));


imagesc3d2(abs(Lpf), N/2, 41, [90,90,90], [0,0.02], 0, ['Lap Ph Full'])
imagesc3d2(abs(Lpx), N/2, 42, [90,90,90], [0,0.02], 0, ['Lap Ph Local'])
imagesc3d2(abs(Lpf-Lpx), N/2, 43, [90,90,90], [0,0.02], 0, ['Lap Ph Back'])
imagesc3d2(abs(Lpu), N/2, 44, [90,90,90], [0,0.08], 0, ['Lap Ph Use'])

imagesc3d2(G, N/2, 1, [90,90,90], [0,1], 0, ['G'])
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
rmse_tvnl = compute_rmse(real(outnl.x.*mask_use),chi);%100 * norm(outnl.x(:).*mask_use(:) - chi(:)) / norm(chi(:));
metrics_nltv = compute_metrics(real(outnl.x.*mask_use),chi);

imagesc3d2(real(outnl.x) .* mask_use - (mask_use==0), N/2, 9, [90,90,90], [-0.12,0.12], 0, ['nlTV RMSE: ', num2str(rmse_tvnl), '  iter : ', num2str(outnl.iter)])


% nlTGV
%params.mu0 = 2*params.mu1;            % second order gradient consistency
%params.alpha0 = 2 * params.alpha1;  % second order gradient L1 penalty

out2nl = nlTGV(params); 
rmse_tgvnl = compute_rmse(real(out2nl.x.*mask_use),chi);%100 * norm(out2nl.x(:).*mask_use(:) - chi(:)) / norm(chi(:));
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
rmse_tvbrw = compute_rmse(real(outbrw.x.*mask_use),chi);%100 * norm(outbrw.x(:).*mask_use(:) - chi(:)) / norm(chi(:));
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

 
 %% Background removal
 
pLBV = LBV(unwrapLaplacian(phase_use2,N,spatial_res),mask_use,N,spatial_res);
 
imagesc3d2(pLBV/10-phase_true, N/2, 50, [90,90,90], [-0.12,0.12], 0, 'pLBV');
imagesc3d2(phase_use-pLBV/10, N/2, 51, [90,90,90], [-0.12,0.12], 0, 'pLBV');


ppdf = [];
ppdf.maxOuterIter = 150;
ppdf.tol_update = 0.01;
ppdf.K = kernel;
ppdf.input = phase_use2;
ppdf.mask = mask_use;
ppdf.weight = magn_use.*mask_use;
ppdf.mu1 = 0.003;                  % gradient consistency
ppdf.alpha1 = 0.0001;%eps;% -0.02;               % gradient L1 penalty
pout = nlPDFb(ppdf);

imagesc3d2(pout.x/10, N/2, 51, [90,90,90], [-0.12,0.12], 0, 'pwPDF');
imagesc3d2(pout.phi/10, N/2, 52, [90,90,90], [-0.12,0.12], 0, 'pwPDF');
imagesc3d2(mask_use.*(phase_full-phase_true-pout.phi/10), N/2, 53, [90,90,90], [-0.12,0.12], 0, 'pwPDF');

ptfi = params;
ptfi.maxOuterIter = 500;
ptfi.tol_update = 0.1;
ptfi.input = phase_use2;
ptfi.background = pout.phi;
ptfi.alpha1 = 10*params.alpha1;
ptfi.mu1 = 10*params.mu1;
ptfi.beta = 10;
ptfi.muh = 100*10;

tfout = pWH_nlTV(ptfi);

imagesc3d2(tfout.x/10, N/2, 54, [90,90,90], [-0.12,0.12], 0, 'x');
imagesc3d2(tfout.phi/10, N/2, 55, [90,90,90], [-0.12,0.12], 0, 'phi');

%muh = 10.^( ((1:17)-9)/2 );%10.^( ((1:25)-15)/4 +2);
mupdf = 10.^( -5+((1:16))/4 );
for sb = 1:16
%for sm = 1:17
%     
ppdf.alpha1 = mupdf(sb);                  % gradient consistency
%ppdf.mu1 = 50*mupdf(sb);                  % gradient consistency
pout = wPDFc(ppdf);

rmsepdfc2(sb) = compute_rmse(real(phase_use-pout.phi).*mask_use,phase_true.*mask_use);
% 
% outss = sswTVc(params); 
% chiss = outss.x;
% 
% rmse3(sb) = compute_rmse(real(chiss.*mask_use),chi_true.*mask_use);
% ssim3(sb) = compute_ssim(real(chiss.*mask_use),chi_true.*mask_use);
clc
display(sb);
end




ppdf.maxOuterIter = 150;
ppdf.tol_update = 0.1;
ppdf.mu2 = 1;                  % gradient consistency
ppdf.alpha1 = 2e-4;% -0.02;               % gradient L1 penalty
ppdf.mu1 = 1e-2;% -0.02;               % gradient L1 penalty
poutl1 = wPDFc(ppdf);

imagesc3d2(poutl1.phi, N/2, 54, [90,90,90], [-0.12,0.12], 0, 'pwPDF');
imagesc3d2(phase_use-poutl1.phi, N/2, 55, [90,90,90], [-0.12,0.12], 0, 'pwPDF');

%ppdf.maxOuterIter = 100;
%ppdf.tol_update = 0.5;
%ppdf.mu1 = 0.001;                  % gradient consistency
ppdf.alpha1 = 0.01;% -0.02;               % gradient L1 penalty
poutl1 = wPDFl1b(ppdf);

imagesc3d2(poutl1.phi, N/2, 56, [90,90,90], [-0.12,0.12], 0, 'pwPDF');
imagesc3d2(phase_use-poutl1.phi, N/2, 57, [90,90,90], [-0.12,0.12], 0, 'pwPDF');
% Common parameters
params2 = [];
params2.K = kernel;
params2.input = pLBV;

%params.mu2 = 1.0;
params2.mu1 = 1e-2;                  % gradient consistency
params2.alpha1 = 2e-4;               % gradient L1 penalty

params2.weight = magn_use;
outw2 = wTV(params2); 
rmse_tvw2 = 100 * norm(outw2.x(:).*mask_use(:) - chi_true(:)) / norm(chi_true(:));
metrics_tvw2 = compute_metrics(real(outw2.x.*mask_use),chi_true);

imagesc3d2(real(outw2.x) .* mask_use - (mask_use==0), N/2, 52, [90,90,90], [-0.12,0.12], 0, ['Weighted TV2 RMSE: ', num2str(rmse_tvw2), '  iter : ', num2str(outw2.iter)])

params2.mask = mask_use;
params2.beta = 50;
params2.muh = 100*50;
outss2 = sswTV(params2); 
chiss2 = outss2.x;
rmse_tvss2 = compute_rmse(real(chiss2.*mask_use),chi_true.*mask_use);
imagesc3d2(real(outss2.x) .* mask_use - (mask_use==0), N/2, 53, [90,90,90], [-0.12,0.12], 0, ['SS TV RMSE: ', num2str(rmse_tvss2), '  iter : ', num2str(outss2.iter)])
imagesc3d2(real(outss2.phi), N/2, 54, [90,90,90], [-0.12,0.12], 0, ['SS Phi'])

[lpPDF hpPDF] = PDF(phase_use, 1.0./(1000*magn_use+eps), mask_use,N,spatial_res, [0 0 1]);


imagesc3d2(phase_use-lpPDF, N/2, 55, [90,90,90], [-0.12,0.12], 0, 'PDF');


params2 = [];
params2.K = kernel;
params2.input = lpPDF;

params2.mu1 = 1e-2;                  % gradient consistency
params2.alpha1 = 2e-4;               % gradient L1 penalty

params2.weight = magn_use;
outw3 = wTV(params2); 
rmse_tvw3 = 100 * norm(outw3.x(:).*mask_use(:) - chi_true(:)) / norm(chi_true(:));
metrics_tvw3 = compute_metrics(real(outw3.x.*mask_use),chi_true);

imagesc3d2(real(outw3.x) .* mask_use - (mask_use==0), N/2, 56, [90,90,90], [-0.12,0.12], 0, ['Weighted TV3 RMSE: ', num2str(rmse_tvw3), '  iter : ', num2str(outw3.iter)])

params2.maxOuterIter = 500;
params2.tol_update = 0.01;
params2.mask = mask_use;
params2.beta = 30;
params2.muh = 100*30;
outss3 = sswTV(params2); 
chiss3 = outss3.x;
rmse_tvss3 = compute_rmse(real(chiss3.*mask_use),chi_true.*mask_use);
imagesc3d2(real(outss3.x) .* mask_use - (mask_use==0), N/2, 57, [90,90,90], [-0.12,0.12], 0, ['SS TV RMSE: ', num2str(rmse_tvss3), '  iter : ', num2str(outss3.iter)])
imagesc3d2(real(outss3.phi), N/2, 58, [90,90,90], [-0.12,0.12], 0, ['SS Phi'])


params3 = params;

params3.mask = mask_use;
params3.maxOuterIter = 100;
params3.tol_update = 0.5;
params3.beta = 100;
params3.muh = 100*100;
outms = sswTV(params3); 
rmse_ms(1) = compute_rmse(real(outms.x.*mask_use),chi_true.*mask_use);
ssim_ms(1) = compute_ssim(real(outms.x.*mask_use),chi_true.*mask_use);
imagesc3d2(real(outms.x) .* mask_use - (mask_use==0), N/2, 59, [90,90,90], [-0.12,0.12], 0, ['MS TV RMSE: ', num2str(rmse_ms(1)), '  iter : ', num2str(outms.iter)])
imagesc3d2(real(outms.phi), N/2, 60, [90,90,90], [-0.12,0.12], 0, ['MS Phi'])

phi(:,:,:,1) = outms.phi;

for c = 2:15
params3.beta = params3.beta/2;
params3.muh = params3.muh/2;
params3.input = params3.input-outms.phi;
outms = sswTV(params3); 

rmse_ms(c) = compute_rmse(real(outms.x.*mask_use),chi_true.*mask_use);
ssim_ms(c) = compute_ssim(real(outms.x.*mask_use),chi_true.*mask_use);
phi(:,:,:,c) = outms.phi;
end

imagesc3d2(real(outms.x) .* mask_use - (mask_use==0), N/2, 59, [90,90,90], [-0.12,0.12], 0, ['MS TV RMSE: ', num2str(rmse_ms(c)), '  iter : ', num2str(outms.iter)])
imagesc3d2(real(outms.phi), N/2, 60, [90,90,90], [-0.12,0.12], 0, ['MS Phi'])


dphi = phi(:,:,:,1);
phi(:,:,:,1) = outms.phi;

for c = 2:9
    dphi = dphi+phi(:,:,:,c);
end


params3 = params;
params3.maxOuterIter = 250;
params3.tol_update = 0.05;
params3.beta = 25;
params3.muh = 100*25;
params3.input = phase_use-pout.phi;
outw3 = sswTV(params3); 
rmse_tvw3 = 100 * norm(outw3.x(:).*mask_use(:) - chi_true(:)) / norm(chi_true(:));

imagesc3d2(real(outw3.x) .* mask_use - (mask_use==0), N/2, 61, [90,90,90], [-0.12,0.12], 0, ['Weighted TV3 RMSE: ', num2str(rmse_tvw3), '  iter : ', num2str(outw3.iter)])





%muh = 10.^( ((1:17)-9)/2 );%10.^( ((1:25)-15)/4 +2);
gamma = 10.^( ((1:15)-24)/4 );
 params.beta = 30;
 params.muh = 100*30;%muh(sm);
for sb = 1:15
%     
 params.gamma = gamma(sb);
outss = sswTVb(params); 
chiss = outss.x;

rmse2g(sb) = compute_rmse(real(chiss.*mask_use),chi_true.*mask_use);
ssim2g(sb) = compute_ssim(real(chiss.*mask_use),chi_true.*mask_use);

clc
display(sb);
end
params.gamma = 2e-3;
outss = sswTVb(params); 
chiss = outss.x;

compute_rmse(real(chiss.*mask_use),chi_true.*mask_use)
compute_ssim(real(chiss.*mask_use),chi_true.*mask_use)
