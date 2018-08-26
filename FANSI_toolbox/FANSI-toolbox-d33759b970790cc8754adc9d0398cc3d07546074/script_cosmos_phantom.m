% This is a sample script to show how to use the functions in this toolbox, 
% and how to set the principal variables.
% This example uses a brain phantom based on a COSMOS reconstruction.
%
% Based on the code by Bilgic Berkin at http://martinos.org/~berkin/software.html
% Last modified by Carlos Milovic in 2017.12.27
%

set(0,'DefaultFigureWindowStyle','docked')
addpath(genpath(pwd))

%% load COSMOS data

load spatial_res;           % voxel size
load msk;                   % brain mask => obtained by eroding the BET mask by 5 voxels (by setting peel=5 in LBV)


load magn;                  % magnitude from transversal orientation
load chi_cosmos;            % COSMOS from 12 orientations (in ppm)

N = size(chi_cosmos);

center = N/2 + 1;

TE = 25e-3;
B0 = 2.8936;
gyro = 2*pi*42.58;

phs_scale = TE * gyro * B0;

imagesc3d2(chi_cosmos, N/2, 1, [90,90,-90], [-0.10,0.14], [], 'COSMOS')


%% Create dipole kernel and susceptibility to field model

kernel = dipole_kernel( N, spatial_res, 0 ); % 0 for the continuous kernel by Salomir and Marques and Bowtell.
                                             % 1 for a discrete kernel
                                             % formulation
                                             % 2 for an integrated Green
                                             % function


% Create raw phase
rphase = ifftn(fftn(chi_cosmos).*kernel);

signal = magn.*exp(1i*phs_scale*(rphase))+0.025*(randn(size(magn))+1i*randn(size(magn))); %SNR = 40
phase_use = angle(signal)/phs_scale;
magn_use = magn;
mask_use = msk;

clear magn msk signal rphase
imagesc3d2(magn_use, N/2, 2, [90,90,-90], [0,1], [], 'Magn')
imagesc3d2(mask_use, N/2, 3, [90,90,-90], [0,1], [], 'Mask')
imagesc3d2(phase_use, N/2, 4, [90,90,-90], [-0.10,0.14], [], 'Local Phase')

%% Reconstruction examples 
% Weighted Linear TV and TGV ADMM reconstructions


% Common parameters
params = [];

params.K = kernel;
params.input = mask_use.*phase_use*phs_scale; % Linear may work with PPM
% We use radians to mantain consistancy with nonlinear parameters


alpha = 8e-3;
params.alpha1 = alpha;
params.mu1 = 10*alpha;

%params.alpha1 = 10^(-(27+10)/10);                 % gradient L1 penalty
%params.mu1 = 10*params.alpha1;                  % gradient consistency

% Linear, no weight
params.weight = ones(N);

% TV
out = wTV(params);
metrics_tv = compute_metrics(real(out.x.*mask_use/phs_scale),chi_cosmos);
imagesc3d2(mask_use.*out.x/phs_scale, N/2, 8, [90,90,-90], [-0.10,0.14], [], 'QSM: TV')

% TGV
out2 = wTGV(params);
metrics_tgv = compute_metrics(real(out2.x.*mask_use/phs_scale),chi_cosmos);
imagesc3d2(mask_use.*out2.x/phs_scale, N/2, 9, [90,90,-90], [-0.10,0.14], [], 'QSM: TGV')


% Linear, magnitude weighted
alpha = 4e-4;
params.alpha1 = alpha;
params.mu1 = 10*alpha;
params.weight = magn_use;

% TV
wout = wTV(params);
metrics_wtv = compute_metrics(real(out.x.*mask_use/phs_scale),chi_cosmos);
imagesc3d2(mask_use.*wout.x/phs_scale, N/2, 12, [90,90,-90], [-0.10,0.14], [], 'QSM: Weighted TV')

% TGV
wout2 = wTGV(params);
metrics_wtgv = compute_metrics(real(out2.x.*mask_use/phs_scale),chi_cosmos);
imagesc3d2(mask_use.*wout2.x/phs_scale, N/2, 13, [90,90,-90], [-0.10,0.14], [], 'QSM: Weighted TGV')


% Nonlinear reconstructions
%params.alpha1 = 1.0*10^(-(28+20)/10);
%params.mu1 = 1.0*10^(-28/10);

outnl = nlTV(params);
chi_nl = outnl.x/phs_scale;
metrics_nltv = compute_metrics(real(chi_nl.*mask_use),chi_cosmos);
imagesc3d2(mask_use.*chi_nl, N/2, 14, [90,90,-90], [-0.10,0.14], [], 'QSM: nTV')

outnl2 = nlTGV(params);
chi_nl2 = outnl2.x/phs_scale;
metrics_nltgv = compute_metrics(real(chi_nl2.*mask_use),chi_cosmos);
imagesc3d2(mask_use.*chi_nl2, N/2, 15, [90,90,-90], [-0.10,0.14], [], 'QSM: nTGV')



%% Optimize for RMSE Example


params = [];

params.K = kernel;
params.input = mask_use.*phase_use*phs_scale; % Linear may work with PPM

for ka = 1:30
    
    alpha = 10^(-1.5-ka/10);
    params.alpha1 = alpha;
    params.mu1 = 10*alpha;                  % gradient consistency
    
    
    params.weight = ones(N);
    out = wTV(params);
    chi = out.x/phs_scale;
    rmse_l2(ka) = compute_rmse(real(chi.*mask_use),chi_cosmos);
    
    params.weight = magn_use;
    out = wTV(params);
    chi = out.x/phs_scale;
    rmse_w(ka) = compute_rmse(real(chi.*mask_use),chi_cosmos);
    
    out = nlTV(params);
    chi = out.x/phs_scale;
    rmse_nl(ka) = compute_rmse(real(chi.*mask_use),chi_cosmos);
end



figure(20);
plot(log(rmse_l2),'g');
hold on;
plot(log(rmse_w),'r');
plot(log(rmse_nl),'b');
hold off;

[kval, ka] = min(rmse_l2);
params.weight = ones(size(phase_use));
alpha = 10^(-1.5-ka/10);
params.alpha1 = alpha;
params.mu1 = 10*alpha;                  % gradient consistency
out = wTV(params);
chi = out.x/phs_scale;

[kval, ka] = min(rmse_w);
params.weight = magn_use;
alpha = 10^(-1.5-ka/10);
params.alpha1 = alpha;
params.mu1 = 10*alpha;                  % gradient consistency
outw = wTV(params);
chiw = outw.x/phs_scale;
outnl = wTV(params);
chinl = outnl.x/phs_scale;


    
   
    