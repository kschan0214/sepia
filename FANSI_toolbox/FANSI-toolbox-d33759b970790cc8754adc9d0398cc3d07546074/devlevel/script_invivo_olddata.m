% This is a sample script to show how to use the functions in this toolbox, 
% and how to set the principal variables.
% This example uses a data set from the Wellcome Trust Centre for
% NeuroImaging, University College London, UK.
%
% Last modified by Carlos Milovic in 2017.06.12
%%%-------------------------------------------------------------------------
%% load data
%%-------------------------------------------------------------------------


set(0,'DefaultFigureWindowStyle','docked')

addpath(genpath(pwd))

load RDF.mat


N = size(Mask);



imagesc3d2(Mask, N/2, 1, [90,90,-90], [0,1], [], 'Mask')


imagesc3d2(RDF, N/2, 2, [90,90,-90], [-1,1], [], 'Input Phase')
imagesc3d2(iMag, N/2, 3, [90,90,-90], [0,0.5], [], 'Magnitude')



%%-------------------------------------------------------------------------
%% create dipole kernel
%%-------------------------------------------------------------------------
spatial_res = voxel_size;

kernel = dipole_kernel_fansi( N, spatial_res, 0, B0_dir );

phs_use = RDF;
magn_use = Mask.*iMag/max(iMag(:));
imagesc3d2(magn_use, N/2, 4, [90,90,-90], [0,1], [], 'rMagnitude')
imagesc3d2(N_std, N/2, 5, [90,90,-90], [0,1], [], 'noise')

%% Phase to PPM scale


TE = 9.09e-3;
B0 = 6.979;
gyro = 2*pi*42.58;

phs_scale = TE * gyro * B0;



%% L-curve example

params = [];

params.input = phs_use;
params.K = kernel;
params.maxOuterIter = 75;

% options = [];
% options.tgv = false;
% options.nonlinear = true;

alpha = 10.^(-(1:25)/5-1);

for ka = 11:20
    
    params.alpha1 = alpha(ka);
    params.mu1 = 25*alpha(ka);    
    
%     
%     params.weight = ones(size(iMag));
%     out = wTV(params); 
%     
%     [ data_cost, reg_cost ] = compute_costs( out.x.*Mask, params.input.*Mask, kernel );    
%     dc(ka) = data_cost;
%     rc(ka) = reg_cost;
%     
%     chi_tv = out.x*delta_TE;
%     rmse_tv2(ka) = compute_rmse( chi_tv.*mask, QSMc.*mask/(pi*phs_scale));
%     ssim_tv(ka) = compute_ssim( chi_tv.*mask, QSMc.*mask/(pi*phs_scale));
%     
    params.weight = magn_use;    
%     out = wTV(params); 
%     
%     [ data_cost, reg_cost ] = compute_costs( out.x.*Mask, params.input.*Mask, kernel );    
%     dcw(ka) = data_cost;
%     rcw(ka) = reg_cost;
    


    out = wTV(params); 
    chi_nltv = out.x*delta_TE;
%     rmse_nltv(ka) = compute_rmse( chi_nltv.*mask, QSMc.*mask/(pi*phs_scale));
    
    rmse_nltv2(ka-10) = compute_rmse( chi_nltv.*mask, QSMc.*mask/(pi*phs_scale));
    ssim_nltv(ka-10) = compute_ssim( chi_nltv.*mask, QSMc.*mask/(pi*phs_scale));
    
%     
%     [ data_cost, reg_cost ] = compute_costs( out.x.*Mask, params.input.*Mask, kernel );    
%     dcnl(ka) = data_cost;
%     rcnl(ka) = reg_cost;
    
%     
%     
%     outf = FANSI( phs_use, magn_use, spatial_res, params.alpha1, params.mu1, 0.0054, options, B0_dir );
%     
%     [ data_cost, reg_cost ] = compute_costs( outf.x.*Mask, phs_use.*Mask, kernel );    
%     dcf(ka) = data_cost;
%     rcf(ka) = reg_cost;
end

% L-Curve analysis
[ Kappa ] = draw_lcurve( alpha(1:25), rc, dc, 11 );
% [ Kappaw ] = draw_lcurve( alpha(1:17), rcw(1:17), dcw(1:17), 12 );
[ Kappanl ] = draw_lcurve( alpha(1:15), rcnl(1:15), dcnl(1:15), 13 );
% [ Kappaf ] = draw_lcurve( alpha(1:17), rcf(1:17), dcf(1:17), 14 );


% Outcome examples
% outf = FANSI( phs_use, magn_use, spatial_res, 0.0015, 100*0.0015, 0.0054, options, B0_dir );%alpha(9)

    params.alpha1 = alpha(6);
    params.mu1 = 25*alpha(6);    
    params.weight = ones(size(iMag));
    out = wTV(params); 
    chi_tv = out.x*delta_TE;

    params.alpha1 = alpha(14);
    params.mu1 = 25*alpha(14);    
    params.weight = magn_use;
    out = nlTV(params); 
    chi_nltv = out.x*delta_TE;

chi_medi = QSMc/(pi*phs_scale);
imagesc3d2(mask.*RDF*delta_TE, N/2, 20, [90,90,-90], [-0.12,0.14], [], 'phase')
imagesc3d2(mask.*chi_tv, N/2, 21, [90,90,-90], [-0.12,0.14], [], 'tv')
imagesc3d2(mask.*chi_nltv, N/2, 22, [90,90,-90], [-0.12,0.14], [], 'nltv')
imagesc3d2(mask.*chi_medi, N/2, 23, [90,90,-90], [-0.12,0.14], [], 'medi')

imagesc3d2(mask.*(chi_nltv-chi_tv), N/2, 24, [90,90,-90], [-0.10,0.10], [], 'dif1')
imagesc3d2(mask.*(chi_nltv-chi_medi), N/2, 25, [90,90,-90], [-0.10,0.10], [], 'dif2')

    params.alpha1 = alpha(6);
    params.mu1 = 25*alpha(6);    
    params.alpha0 = 2*alpha(6);
    params.mu0 = 2*25*alpha(6);    
    params.weight = ones(size(iMag));
    out = wTGV(params); 
    chi_tgv = out.x*delta_TE;

    params.alpha1 = alpha(14);
    params.mu1 = 25*alpha(14);    
    params.alpha0 = 2*alpha(14);
    params.mu0 = 2*25*alpha(14);  
    params.weight = magn_use;
    out = nlTGV(params); 
    chi_nltgv = out.x*delta_TE;
    
    
imagesc3d2(mask.*RDF*delta_TE, N/2, 20, [90,90,-90], [-0.12,0.14], [], 'phase')
imagesc3d2(mask.*chi_tgv, N/2, 21, [90,90,-90], [-0.12,0.14], [], 'tgv')
imagesc3d2(mask.*chi_nltgv, N/2, 22, [90,90,-90], [-0.12,0.14], [], 'nltgv')
imagesc3d2(mask.*chi_medi, N/2, 23, [90,90,-90], [-0.12,0.14], [], 'medi')

imagesc3d2(mask.*(chi_nltgv-chi_tgv), N/2, 24, [90,90,-90], [-0.10,0.10], [], 'dif1')
imagesc3d2(mask.*(chi_nltgv-chi_medi), N/2, 25, [90,90,-90], [-0.10,0.10], [], 'dif2')
imagesc3d2(mask.*(chi_nltgv-chi_nltv), N/2, 26, [90,90,-90], [-0.10,0.10], [], 'dif3')