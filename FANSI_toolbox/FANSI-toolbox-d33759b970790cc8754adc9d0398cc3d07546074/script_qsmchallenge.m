% This is a sample script to show how to use the functions in this toolbox, 
% and how to set the principal variables.
% This example uses a data from the QSM reconstruction challenge 2016.
% 4th International Workshop on Phase Contrast and QSM:
% http://qsm.rocks
%
% Based on the code by Bilgic Berkin at http://martinos.org/~berkin/software.html
% Last modified by Carlos Milovic in 2017.12.27
%%%-------------------------------------------------------------------------
%% load data
%%-------------------------------------------------------------------------


set(0,'DefaultFigureWindowStyle','docked')

addpath(genpath(pwd))


load phs_tissue;            % tissue phase from transversal orientation (in ppm, normalized by gyro*TE*B0)
load spatial_res;           % voxel size
load msk;                   % brain mask => obtained by eroding the BET mask by 5 voxels (by setting peel=5 in LBV)

load magn;                  % magnitude from transversal orientation

N = size(msk);



imagesc3d2(msk, N/2, 1, [90,90,-90], [0,1], [], 'Mask')

imagesc3d2(phs_tissue, N/2, 2, [90,90,-90], [-0.05,0.05], [], 'Input Phase')
imagesc3d2(magn, N/2, 3, [90,90,-90], [0,0.5], [], 'Magnitude')


%%-------------------------------------------------------------------------
%% create dipole kernel
%%-------------------------------------------------------------------------

kernel = dipole_kernel( N, spatial_res, 0 );

%%-------------------------------------------------------------------------
%% TKD recon
%%-------------------------------------------------------------------------

thre_tkd = 0.19;      % TKD threshold parameter
chi_tkd = tkd( phs_tissue, msk, kernel, thre_tkd, N );

imagesc3d2(chi_tkd, N/2, 11, [90,90,-90], [-0.10,0.14], [], 'TKD')


%%-------------------------------------------------------------------------
%% closed-form L2 recon
%%-------------------------------------------------------------------------

l2beta = 9e-2;    % regularization parameter
chi_L2 = chiL2( phs_tissue, msk, kernel, l2beta, N );

imagesc3d2(chi_L2, N/2, 12, [90,90,-90], [-0.10,0.14], [], 'CF L2')

%% Linear vs Weighted Linear TV

num_iter = 50;
tol_update = 1;

params = [];

params.mu1 = 1e-4;                  % gradient consistency
params.alpha1 = 2e-4;               % gradient L1 penalty

params.maxOuterIter = num_iter;
params.tol_update = tol_update;

params.K = kernel;
params.input = phs_tissue;
 

out = wTV(params); 

magn_use = magn .* msk;
magn_use = sum(msk(:))*magn_use / sum(abs(magn_use(:)));

params.weight = magn_use;
outw = wTV(params); 
chiw = outw.x;

imagesc3d2(out.x.*msk, N/2, 21, [90,90,-90], [-0.10,0.14], [], 'TV')
imagesc3d2(chiw.*msk, N/2, 22, [90,90,-90], [-0.10,0.14], [], 'wTV')
    
    
 %% Nonlinear TV
    
TE = 25e-3;
B0 = 2.8936;
gyro = 2*pi*42.58;
phs_scale = TE * gyro * B0;
 
params = [];
params.input = phs_tissue * phs_scale;
magn_use = magn .* msk;
magn_use = magn_use / max(abs(magn_use(:)));
params.weight = magn_use;
params.K = kernel;

params.alpha1 = 4e-4;               % gradient L1 penalty
params.mu1 = 1e-2;                  % gradient consistency

outw2 = wTV(params);
chiw2 = outw2.x/phs_scale;

outnl = nlTV(params); 
chinl = outnl.x/phs_scale;

imagesc3d2(chiw2.*msk, N/2, 23, [90,90,-90], [-0.1,0.14], [], 'wTV2')
imagesc3d2(chinl.*msk, N/2, 24, [90,90,-90], [-0.1,0.14], [], 'nlTV')


%% L-curve example
load chi_cosmos
load chi_sti

params = [];
params.input = phs_tissue * phs_scale;
params.K = kernel;
magn_use = magn .* msk;
magn_use = magn_use / max(abs(magn_use(:)));%norm(abs(magn_use(:)));

for ka = 1:30
    
    alpha(ka) = 10^(-1.5-ka/10);
    params.alpha1 = alpha(ka);
    params.mu1 = 25*alpha(ka);                  % gradient consistency
    
    
    params.weight = ones(size(magn));
    out = wTV(params); 
    
    [ data_cost, reg_cost ] = compute_costs( out.x.*msk, params.input.*msk, kernel );    
    dc(ka) = data_cost;
    rc(ka) = reg_cost;
    
    params.weight = magn_use;    
    out = wTV(params); 
    
    [ data_cost, reg_cost ] = compute_costs( out.x.*msk, params.input.*msk, kernel );    
    dcw(ka) = data_cost;
    rcw(ka) = reg_cost;
    
    % Error evaluation with respect to different ground-truths
    r_cosmos(ka) = compute_rmse(out.x.*msk/phs_scale,chi_cosmos);
    r_x33(ka) = compute_rmse(out.x.*msk/phs_scale,chi_33);
    r_sti(ka) = compute_rmse(out.x.*msk/phs_scale,chi_sti);
    r_stil2(ka) = compute_rmse(out.x.*msk/phs_scale,chi_sti_l2);
    
    s_cosmos(ka) = compute_ssim(out.x.*msk/phs_scale,chi_cosmos);
    s_x33(ka) = compute_ssim(out.x.*msk/phs_scale,chi_33);
    s_sti(ka) = compute_ssim(out.x.*msk/phs_scale,chi_sti);
    s_stil2(ka) = compute_ssim(out.x.*msk/phs_scale,chi_sti_l2); % only for the wTV result, for simplicity
    
    out = nlTV(params); 
    
    [ data_cost, reg_cost ] = compute_costs( out.x.*msk, params.input.*msk, kernel );    
    dcnl(ka) = data_cost;
    rcnl(ka) = reg_cost;
end

[ Kappa ] = draw_lcurve( alpha, rc, dc, 31 );
[ Kappaw ] = draw_lcurve( alpha, rcw, dcw, 32 );
[ Kappanl ] = draw_lcurve( alpha, rcnl, dcnl, 33 );


figure(28)
semilogx(alpha,r_cosmos,'b');
hold on;
semilogx(alpha,r_sti,'r');
semilogx(alpha,r_stil2,'g');
semilogx(alpha,r_x33,'y');
hold off;

figure(51)
semilogx(alpha,s_cosmos,'b');
hold on;
semilogx(alpha,s_sti,'r');
semilogx(alpha,s_stil2,'g');
semilogx(alpha,s_x33,'y');
hold off;

options = [];
options.tgv = false;
options.nonlinear = false;

for ka = 1:30
    
    alpha(ka) = 10^(-1.5-ka/10);
    alpha1 = alpha(ka);
    mu1 = 25*alpha1;                  % gradient consistency
    
    outf = FANSI( phs_tissue * phs_scale, magn_use, spatial_res, alpha1, mu1, 1/40, options, [0 0 1] );
    
    [ data_cost, reg_cost ] = compute_costs( outf.x.*msk, phs_tissue * phs_scale.*msk, kernel );    
    dcf(ka) = data_cost;
    rcf(ka) = reg_cost;
end

[ Kappaf ] = draw_lcurve( alpha, rcf, dcf, 34 );

outf = FANSI( phs_tissue * phs_scale, magn_use, spatial_res, alpha(17), 25*alpha(17), 1/40 );
imagesc3d2(outf.x.*msk/phs_scale, N/2, 25, [90,90,-90], [-0.1,0.14], [], 'x1f')

