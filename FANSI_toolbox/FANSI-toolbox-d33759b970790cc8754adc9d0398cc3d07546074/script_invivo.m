% This is a sample script to show how to use the functions in this toolbox, 
% and how to set the principal variables.
% This example uses a data set from the Wellcome Trust Centre for
% NeuroImaging, University College London, UK.
%
% Last modified by Carlos Milovic in 2017.12.27
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

TE = 20e-3;
B0 = 3;
gyro = 2*pi*42.58;

phs_scale = TE * gyro * B0;



%% L-curve example

params = [];

params.input = phs_use;
params.K = kernel;
params.maxOuterIter = 50;

options = [];
options.tgv = false;
options.nonlinear = true;
alpha = 10.^(-(1:17)/5-1);
for ka = 1:17
    
    params.alpha1 = alpha(ka);
    params.mu1 = 100*alpha(ka);    
    
    
    params.weight = ones(size(iMag));
    out = wTV(params); 
    
    [ data_cost, reg_cost ] = compute_costs( out.x.*Mask, params.input.*Mask, kernel );    
    dc(ka) = data_cost;
    rc(ka) = reg_cost;
    
    
    params.weight = magn_use;    
    out = wTV(params); 
    
    [ data_cost, reg_cost ] = compute_costs( out.x.*Mask, params.input.*Mask, kernel );    
    dcw(ka) = data_cost;
    rcw(ka) = reg_cost;
%     
%     out = nlTV(params); 
%     
%     [ data_cost, reg_cost ] = compute_costs( out.x.*Mask, params.input.*Mask, kernel );    
%     dcnl(ka) = data_cost;
%     rcnl(ka) = reg_cost;
    
    
    
    outf = FANSI( phs_use, magn_use, spatial_res, params.alpha1, params.mu1, 0.0054, options, B0_dir );
    
    [ data_cost, reg_cost ] = compute_costs( outf.x.*Mask, phs_use.*Mask, kernel );    
    dcf(ka) = data_cost;
    rcf(ka) = reg_cost;
end

% L-Curve analysis
[ Kappa ] = draw_lcurve( alpha, rc, dc, 11 );
[ Kappaw ] = draw_lcurve( alpha(1:17), rcw(1:17), dcw(1:17), 12 );
%[ Kappanl ] = draw_lcurve( alpha(1:17), rcnl(1:17), dcnl(1:17), 13 );
[ Kappaf ] = draw_lcurve( alpha(1:17), rcf(1:17), dcf(1:17), 14 );


% Outcome examples
outf = FANSI( phs_use, magn_use, spatial_res, 0.0015, 100*0.0015, 0.0054, options, B0_dir );%alpha(9)


imagesc3d2(outf.x, N/2, 21, [90,90,-90], [-1,1], [], 'Xfansi optimal')

outfc1 = FANSI( phs_use, magn_use, spatial_res, alpha(16), 100*alpha(16), 0.0054, options, B0_dir );
outfc2 = FANSI( phs_use, magn_use, spatial_res, alpha(17), 100*alpha(17), 0.0054, options, B0_dir );


imagesc3d2(outfc1.x, N/2, 22, [90,90,-90], [-1,1], [], 'Xfansi c1')
imagesc3d2(outfc2.x, N/2, 23, [90,90,-90], [-1,1], [], 'Xfansi c2')

outfc3 = FANSI( phs_use, magn_use, spatial_res, alpha(13), 100*alpha(13), 0.0054, options, B0_dir );
imagesc3d2(outfc3.x, N/2, 24, [90,90,-90], [-1,1], [], 'Xfansi c3')
outfc4 = FANSI( phs_use, magn_use, spatial_res, alpha(11), 100*alpha(11), 0.0054, options, B0_dir );
imagesc3d2(outfc4.x, N/2, 26, [90,90,-90], [-1,1], [], 'Xfansi c4')