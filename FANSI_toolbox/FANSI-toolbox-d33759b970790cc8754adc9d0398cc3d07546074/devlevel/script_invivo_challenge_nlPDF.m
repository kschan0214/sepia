% This is a sample script to show how to use the functions in this toolbox, 
% and how to set the principal variables.
% This example uses a brain phantom based on a COSMOS reconstruction.
%
% Based on the code by Bilgic Berkin at http://martinos.org/~berkin/software.html
% Last modified by Carlos Milovic in 2017.06.12
%

set(0,'DefaultFigureWindowStyle','docked')

addpath(genpath(pwd))

%% load COSMOS data
load('test_cosmos_nlPDF.mat')

TE = 25e-3;
B0 = 2.8936;
gyro = 2*pi*42.58;

phs_scale = TE * gyro * B0;

%% Create dipole kernel and susceptibility to field model

kernel = dipole_kernel( N, spatial_res, 0 ); % 0 for the continuous kernel by Salomir and Marques and Bowtell.
                                             % 1 for a discrete kernel
                                             % formulation
                                             % 2 for an integrated Green
                                             % function

kernelg = dipole_kernel( N, spatial_res, 2 );
% Create raw phase
phase_use = angle(Img);
magn_use = magn/max(magn(:));
magn = Magn_orig/0.033;
mask_use = msk;
mask_bet = msk+msk_deleted;
mask_full = max(real(mp_rage>0.02),mask_bet);%msk+msk_deleted+msk_tissue;

se = strel('sphere',1);
mask_full=imdilate(mask_full,se);
mask_full=imerode(mask_full,se);

imagesc3d2(phase_use, N/2, 10, [90,90,-90], [-3.14,3.14], [], 'Phase')
imagesc3d2(angle(Img), N/2, 11, [90,90,-90], [-3.14,3.14], [], 'Phase single acq')
imagesc3d2(angle(Img.*exp(1i*pi*0.5))-angle(Img), N/2, 12, [90,90,-90], [-3.14,3.14], [], 'Phase single acq')
imagesc3d2(mask_bet, N/2, 3, [90,90,-90], [0,1], [], 'Mask')
imagesc3d2(mask_full, N/2, 4, [90,90,-90], [0,1], [], 'Mask')
imagesc3d2(mp_rage/max(mp_rage(:)), N/2, 5, [90,90,-90], [0,0.04], [], 'Magn')
imagesc3d2(magn, N/2, 6, [90,90,-90], [0,1], [], 'Magn')



%% Unwrapping results
phase_unwrap = unwrap(phase_use,spatial_res);
phase_unwrapL = unwrapLaplacian(phase_use,N,spatial_res);

phase_unwrapL2 = unwrapLaplacian(angle(exp(1i*(phase_use-phase_unwrapL))),N,spatial_res);
phase_unwrapL3 = unwrapLaplacian(angle(exp(1i*(phase_use-phase_unwrapL-phase_unwrapL2))),N,spatial_res);
phase_unwrapL4 = unwrapLaplacian(angle(exp(1i*(phase_use-phase_unwrapL-phase_unwrapL2-phase_unwrapL3))),N,spatial_res);
phase_unwrapL5 = unwrapLaplacian(angle(exp(1i*(phase_use-phase_unwrapL-phase_unwrapL2-phase_unwrapL3-phase_unwrapL4))),N,spatial_res);

imagesc3d2(phase_unwrapL, N/2, 21, [90,90,-90], [-6,6], [], 'PuwL')
imagesc3d2(phase_unwrap, N/2, 22, [90,90,-90], [-6,6], [], 'Puw')


imagesc3d2(phase_unwrapL2.*mask_bet, N/2, 23, [90,90,-90], [-6,6], [], 'PuwL2')
imagesc3d2(phase_unwrapL3.*mask_bet, N/2, 24, [90,90,-90], [-6,6], [], 'PuwL3')
imagesc3d2(phase_unwrapL4.*mask_bet, N/2, 25, [90,90,-90], [-6,6], [], 'PuwL4')
imagesc3d2(phase_unwrapL5.*mask_bet, N/2, 25, [90,90,-90], [-6,6], [], 'PuwL5')


imagesc3d2(angle(exp(1i*(phase_use-phase_unwrapL-phase_unwrapL2-phase_unwrapL3-phase_unwrapL4))).*mask_bet, N/2, 30, [90,90,-90], [-6,6], [], 'res phase')
imagesc3d2(angle(exp(1i*(phase_use-phase_unwrapL-phase_unwrapL2-phase_unwrapL3-phase_unwrapL4-phase_unwrapL5))).*mask_bet, N/2, 31, [90,90,-90], [-6,6], [], 'res phase')


delta = phase_use;
background = zeros(N);

for it = 1:1000
puL = unwrapLaplacian(delta,N,spatial_res);
background = background+puL;
delta = angle(exp(1i*(phase_use-background)));
if max(delta(:))<pi/2 && min(delta(:)) > -pi/2
    break
end
end
imagesc3d2(delta.*mask_bet, N/2, 32, [90,90,-90], [-6,6], [], 'res phase')
imagesc3d2(delta.*mask_bet, N/2, 33, [90,90,-90], [-6,6], [], 'res phase')

phase_Munwrap = unwrap(phase_use.*mask_bet,spatial_res);
phase_MunwrapL = unwrapLaplacian(phase_use.*mask_bet,N,spatial_res);

imagesc3d2(phase_MunwrapL.*mask_bet, N/2, 23, [90,90,-90], [-6,6], [], 'PuwL')
imagesc3d2(phase_Munwrap.*mask_bet, N/2, 24, [90,90,-90], [-6,6], [], 'Puw')

% Unwrapped filtering
load ph_sti;
lphase = phase_sti*phs_scale;
imagesc3d2(nlocalp, N/2, 29, [90,90,-90], [-1.12,1.12], 0, 'local phase');
imagesc3d2(lphase, N/2, 30, [90,90,-90], [-1.12,1.12], 0, 'local phase');

pLBV = LBV(phase_unwrapL,mask_use,N,spatial_res);
imagesc3d2(pLBV, N/2, 31, [90,90,-90], [-1.12,1.12], 0, 'pLBV');
imagesc3d2(pLBV-nlocalp, N/2, 32, [90,90,-90], [-1.12,1.12], 0, 'pLBV');

[lpPDF hpPDF] = PDF(phase_unwrapL, 1./(40*magn_use+eps), mask_use,N,spatial_res, [0 0 1]);
imagesc3d2(lpPDF, N/2, 33, [90,90,-90], [-1.12,1.12], 0, 'PDF');
imagesc3d2(lpPDF-nlocalp, N/2, 34, [90,90,-90], [-1.12,1.12], 0, 'PDF');


pLBVbet = LBV(phase_unwrapL,mask_bet,N,spatial_res);
imagesc3d2(pLBVbet, N/2, 35, [90,90,-90], [-1.12,1.12], 0, 'pLBV');
imagesc3d2(pLBVbet-nlocalp, N/2, 36, [90,90,-90], [-1.12,1.12], 0, 'pLBV');

[lpPDFbet hpPDFbet] = PDF(phase_unwrapL, 1./(40*magn_use+eps), mask_bet,N,spatial_res, [0 0 1]);
imagesc3d2(lpPDFbet, N/2, 37, [90,90,-90], [-1.12,1.12], 0, 'PDF');
imagesc3d2(lpPDFbet-nlocalp, N/2, 38, [90,90,-90], [-1.12,1.12], 0, 'PDF');


pLBVi = LBV(lapunwrap,mask_use,N,spatial_res);
imagesc3d2(pLBVi, N/2, 31, [90,90,-90], [-1.12,1.12], 0, 'pLBV');
imagesc3d2(pLBVi-nlocalp, N/2, 32, [90,90,-90], [-1.12,1.12], 0, 'pLBV');

[lpPDFi hpPDFi] = PDF(lapunwrap, 1./(40*magn_use+eps), mask_use,N,spatial_res, [0 0 1]);
imagesc3d2(lpPDFi, N/2, 33, [90,90,-90], [-1.12,1.12], 0, 'PDF');
imagesc3d2(lpPDFi-nlocalp, N/2, 34, [90,90,-90], [-1.12,1.12], 0, 'PDF');
%% ADMM PDF


% Linear
ppdf = [];
ppdf.maxOuterIter = 150;
ppdf.tol_update = 0.01;
ppdf.K = kernel;
ppdf.input = phase_unwrapL;
ppdf.mask = mask_use;
ppdf.weight = magn_use.*mask_use;
ppdf.mu1 = 0.003;                  % gradient consistency
ppdf.alpha1 = eps;%0.0001;              % gradient L1 penalty
pout = wPDFb(ppdf);

imagesc3d2((phase_unwrapL-pout.phi).*mask_use, N/2, 35, [90,90,-90], [-1.12,1.12], 0, 'wPDF');

ppdf.input = phase_unwrap;
pout2 = wPDFb(ppdf);
imagesc3d2((phase_unwrap-pout2.phi).*mask_use, N/2, 36, [90,90,-90], [-1.12,1.12], 0, 'wPDF');


mupdf = 10.^( -5+((1:16))/4 );
for sb = 1:16
display(sb);
ppdf.mu1 = mupdf(sb); 
out = wPDFb(ppdf);

rmsepdfL(sb) = compute_rmse(real(out.local).*mask_use,lphase.*mask_use);
rmsepdfL2(sb) = compute_rmse(real(out.local).*mask_use,(rphase-lphase).*mask_use);
rmsepdfB(sb) = compute_rmse(real(out.phi).*mask_use,bphase.*mask_use);
clc
end

ppdf.mu1 = mupdf(9); 
out1 = wPDFb(ppdf);
ppdf.mu1 = mupdf(5); 
ou2 = wPDFb(ppdf);
ppdf.mu1 = mupdf(1); 
ou3 = wPDFb(ppdf);
ppdf.mu1 = mupdf(16); 
ou4 = wPDFb(ppdf);

imagesc3d2(out1.local.*mask_use, N/2, 41, [90,90,-90], [-1.12,1.12], 0, 'wPDF');
imagesc3d2(ou2.local.*mask_use, N/2, 42, [90,90,-90], [-1.12,1.12], 0, 'wPDF');
imagesc3d2(ou3.local.*mask_use, N/2, 43, [90,90,-90], [-1.12,1.12], 0, 'wPDF');
imagesc3d2(ou4.local.*mask_use, N/2, 44, [90,90,-90], [-1.12,1.12], 0, 'wPDF');

ppdf = [];
ppdf.maxOuterIter = 1500;
ppdf.tol_update = 0.001;
ppdf.K = kernel;
ppdf.input = phase_unwrapL;
ppdf.mask = mask_use;
ppdf.weight = magn_use.*mask_use;
ppdf.mu1 = 0.0018;                  % gradient consistency
ppdf.alpha1 = eps;%0.0001;              % gradient L1 penalty
ou5 = wPDFb(ppdf);
imagesc3d2(ou5.local.*mask_use, N/2, 45, [90,90,-90], [-1.12,1.12], 0, 'wPDF');
imagesc3d2(ou5.local.*mask_use-nlocalp, N/2, 45, [90,90,-90], [-1.12,1.12], 0, 'wPDF');
imagesc3d2(ou5.local.*mask_use-nlocalp, N/2, 81, [90,90,-90], [-0.12,0.12], 0, 'wPDF');

ppdf.input = phase_unwrap;
ou6 = wPDFb(ppdf);
imagesc3d2(ou6.local.*mask_use, N/2, 46, [90,90,-90], [-1.12,1.12], 0, 'wPDF');


% Nonlinear
imagesc3d2((bphase-3.7479).*mask_use, N/2, 50, [90,90,-90], [-6,6], 0, 'nlPDF');

ppdf = [];
ppdf.maxOuterIter = 150;
ppdf.tol_update = 0.01;
ppdf.K = kernel;
ppdf.input = phase_use;
ppdf.mask = mask_use;
ppdf.weight = magn_use.*mask_use;
ppdf.mu1 = 0.0018;                  % gradient consistency
ppdf.alpha1 = eps;%0.0001;              % gradient L1 penalty
nout = nlPDFb(ppdf);

imagesc3d2(nout.local.*mask_use, N/2, 51, [90,90,-90], [-1.12,1.12], 0, 'nlPDF');

imagesc3d2(nout.phi.*mask_use, N/2, 52, [90,90,-90], [-6,6], 0, 'nlPDF');

ppdf = [];
ppdf.maxOuterIter = 1500;
ppdf.tol_update = 0.01;
ppdf.K = kernel;
ppdf.input = phase_use;
ppdf.mask = mask_use;
ppdf.weight = magn_use.*mask_use;
ppdf.mu1 = 0.0018;                  % gradient consistency
ppdf.alpha1 = eps;%0.0001;              % gradient L1 penalty
nout2 = nlPDFb(ppdf);

imagesc3d2(nout2.local.*mask_use, N/2, 53, [90,90,-90], [-1.12,1.12], 0, 'nlPDF');

imagesc3d2(nout2.phi.*mask_use, N/2, 54, [90,90,-90], [-6,6], 0, 'nlPDF');

ppdf = [];
ppdf.maxOuterIter = 150;
ppdf.tol_update = 0.01;
ppdf.K = kernel;
ppdf.input = phase_use;
ppdf.mask = mask_use;
ppdf.weight = magn_use.*mask_use;
ppdf.mu1 = 0.0018;                  % gradient consistency
ppdf.alpha1 = eps;%0.0001;              % gradient L1 penalty
nout3 = nlPDFp(ppdf);

imagesc3d2(nout3.local.*mask_use, N/2, 55, [90,90,-90], [-1.12,1.12], 0, 'nlPDF');

imagesc3d2(nout3.phi.*mask_use, N/2, 56, [90,90,-90], [-6,6], 0, 'nlPDF');

ppdf = [];
ppdf.maxOuterIter = 1500;
ppdf.tol_update = 0.001;
ppdf.K = kernel;
ppdf.input = phase_use;
ppdf.mask = mask_use;
ppdf.weight = magn_use.*mask_use;
ppdf.mu1 = 0.0018;                  % gradient consistency
ppdf.alpha1 = eps;%0.0001;              % gradient L1 penalty
nout4 = nlPDFp(ppdf);
imagesc3d2(nout4.local.*mask_use, N/2, 57, [90,90,-90], [-1.12,1.12], 0, 'nlPDF');

imagesc3d2(nout4.phi.*mask_use, N/2, 58, [90,90,-90], [-6,6], 0, 'nlPDF');


ppdf = [];
ppdf.maxOuterIter = 130;
ppdf.tol_update = 0.01;
ppdf.K = kernelg;
ppdf.input = phase_use;
ppdf.mask = mask_use;
ppdf.outermask = mask_full;
ppdf.weight = magn.*mask_full;
ppdf.mu1 = 0.0018;                  % gradient consistency
ppdf.alpha1 = eps;%0.0001;              % gradient L1 penalty
nout5 = nlPDFp(ppdf);

imagesc3d2(nout5t.local.*mask_use, N/2, 71, [90,90,-90], [-1.12,1.12], 0, 'nlPDF');
imagesc3d2(nout5.local.*mask_use-nlocalp, N/2, 72, [90,90,-90], [-1.12,1.12], 0, 'nlPDF');


[kunwrap, lapunwrap, delta] = iterlap_unwrap(phase_use, mask_use, magn.*mask_full, 250 );
imagesc3d2(kunwrap, N/2, 81, [90,90,-90], [-40,40], 0, 'kuw');
imagesc3d2(lapunwrap, N/2, 82, [90,90,-90], [-40,40], 0, 'luw');
imagesc3d2(delta, N/2, 83, [90,90,-90], [-6,6], 0, 'delta');
ppdf.ph_unwrap = kapunwrap;

phi_e = backmodel( mask_full, 0.0, 2 );
ppdf.airmodel = phi_e*phs_scale;
imagesc3d2(ppdf.airmodel.*mask_use, N/2, 84, [90,90,-90], [-40,40], 0, 'air');

[ background, kappa ] = fitmodels( kunwrap, mask_use, magn, phi_e*phs_scale );

phi_e2 = backmodel( ones(N), 0.0, 2 );
[ background2, kappa2 ] = fitmodels( kunwrap, mask_use, magn, phi_e2 );
ppdf.maxOuterIter = 1500;
ppdf.tol_update = 0.001;
ppdf.input = phase_unwrapL;%kunwrap-background2;
ppdf.mu1 = 0.0018;                  % gradient consistency
ppdf.alpha1 = eps;%0.0001;              % gradient L1 penalty
pout = wPDFb(ppdf);
imagesc3d2(pout.phi.*mask_use, N/2, 72, [90,90,-90], [-6,6], 0, 'wPDF');
imagesc3d2((phase_unwrapL-pout.phi-nlocalp).*mask_use, N/2, 73, [90,90,-90], [-6,6], 0, 'wPDF');
ppdf.maxOuterIter = 1500;
ppdf.tol_update = 0.001;
ppdf.input = lapunwrap;%kunwrap-background2;
ppdf.mu1 = 0.0018;                  % gradient consistency
ppdf.alpha1 = eps;%0.0001;              % gradient L1 penalty
pout = wPDFb(ppdf);
imagesc3d2(pout.phi.*mask_use, N/2, 74, [90,90,-90], [-6,6], 0, 'wPDF');
imagesc3d2((lapunwrap-pout.phi-nlocalp).*mask_use, N/2, 75, [90,90,-90], [-1.12,1.12], 0, 'wPDF');

ppdf.maxOuterIter = 150;
ppdf.input = phase_use;
ppdf.mask = mask_use;
ppdf.outermask = mask_full;
ppdf.weight = magn.*mask_full;
ppdf.airmodel = background+pout.phi;
nout5t = nlPDFp3(ppdf);

imagesc3d2(nout5t.phi.*mask_use, N/2, 72, [90,90,-90], [-6,6], 0, 'nlPDF');

ppdf = [];
ppdf.maxOuterIter = 1500;
ppdf.tol_update = 0.001;
ppdf.K = kernel;
ppdf.input = phase_use;
ppdf.mask = mask_use;
ppdf.weight = magn_use.*mask_use;
ppdf.mu1 = 0.0018;                  % gradient consistency
ppdf.alpha1 = eps;%0.0001;              % gradient L1 penalty
nout6 = nlPDFp(ppdf);
imagesc3d2(nout6.local.*mask_use, N/2, 73, [90,90,-90], [-1.12,1.12], 0, 'nlPDF');
imagesc3d2(nout6.local.*mask_use-nlocalp, N/2, 73, [90,90,-90], [-1.12,1.12], 0, 'nlPDF');
imagesc3d2(nout6.local.*mask_use-nlocalp, N/2, 82, [90,90,-90], [-0.12,0.12], 0, 'nlPDF');

imagesc3d2(nout6.phi.*mask_use, N/2, 74, [90,90,-90], [-6,6], 0, 'nlPDF');

ppdf = [];
ppdf.maxOuterIter = 7500;
ppdf.tol_update = 0.001;
ppdf.K = kernel;
ppdf.input = phase_use;
ppdf.mask = mask_use;
ppdf.weight = magn_use.*mask_use;
ppdf.mu1 = 0.0018;                  % gradient consistency
ppdf.alpha1 = eps;%0.0001;              % gradient L1 penalty
nout7 = nlPDFp(ppdf);
imagesc3d2(nout7.local.*mask_use, N/2, 75, [90,90,-90], [-1.12,1.12], 0, 'nlPDF');

imagesc3d2(nout7.phi.*mask_use, N/2, 76, [90,90,-90], [-6,6], 0, 'nlPDF');


% Sensibility test

metrics_nlpdf = compute_metrics(real(nout6.local.*mask_use),nlocalp);

ppdf = [];
ppdf.maxOuterIter = 1500;
ppdf.tol_update = 0.001;
ppdf.K = kernel;
ppdf.input = phase_use;
ppdf.mask = mask_use;
ppdf.weight = magn_use.*mask_use;
ppdf.mu1 = 0.0018*10;                  % gradient consistency
ppdf.alpha1 = eps;%0.0001;              % gradient L1 penalty
nout6t = nlPDFp(ppdf);

metrics_nlpdf_up = compute_metrics(real(nout6t.local.*mask_use),nlocalp);

ppdf.mu1 = 0.0018/10;                  % gradient consistency
nout6t = nlPDFp(ppdf);
metrics_nlpdf_dw = compute_metrics(real(nout6t.local.*mask_use),nlocalp);

%% Reconstruction examples 
% Weighted Linear TV and TGV ADMM reconstructions


% Common parameters
params = [];

params.K = kernel;
params.input = nout6.local; % Linear may work with PPM
% We use radians to mantain consistancy with nonlinear parameters


    alpha = 10^(-24/10-1);
    params.alpha1 = alpha;
    params.mu1 = 25*alpha;           
    
%params.alpha1 = 10^(-(27+10)/10);                 % gradient L1 penalty
%params.mu1 = 100*params.alpha1;                  % gradient consistency

% TV
params.weight = magn_use.*mask_use;



params.input = pLBV; % Linear may work with PPM
out_lbv = wTV(params); 
metrics_tv_lbv = compute_metrics(real(out_lbv.x.*mask_use/phs_scale),chi_cosmos);

params.input = lpPDF; % Linear may work with PPM
out_pdf = wTV(params); 
metrics_tv_pdf = compute_metrics(real(out_pdf.x.*mask_use/phs_scale),chi_cosmos);

params.input = ou5.local; % Linear may work with PPM
out_wpdf = wTV(params); 
metrics_tv_wpdf = compute_metrics(real(out_wpdf.x.*mask_use/phs_scale),chi_cosmos);

params.input = ou6.local; % Linear may work with PPM
out_wpdfi = wTV(params); 
metrics_tv_wpdfi = compute_metrics(real(out_wpdfi.x.*mask_use/phs_scale),chi_cosmos);

params.input = nout2.local; % Linear may work with PPM
out_nlpdf = wTV(params); 
metrics_tv_nlpdf = compute_metrics(real(out_nlpdf.x.*mask_use/phs_scale),chi_cosmos);

params.input = nout6.local; % Linear may work with PPM
out_nlpdfp = wTV(params); 
metrics_tv_nlpdfp = compute_metrics(real(out_nlpdfp.x.*mask_use/phs_scale),chi_cosmos);

imagesc3d2(chi_cosmos, N/2, 101, [90,90,-90], [-0.10,0.14], [], 'COSMOS')
imagesc3d2(out_lbv.x/phs_scale, N/2, 102, [90,90,-90], [-0.10,0.14], [], 'QSM: wtv-lbv')
imagesc3d2(out_pdf.x/phs_scale, N/2, 103, [90,90,-90], [-0.10,0.14], [], 'QSM: wtv-pdf')
imagesc3d2(out_wpdf.x/phs_scale, N/2, 104, [90,90,-90], [-0.10,0.14], [], 'QSM: wtv-wpdf')
imagesc3d2(out_wpdfi.x/phs_scale, N/2, 105, [90,90,-90], [-0.10,0.14], [], 'QSM: wtv-wpdfi')
%imagesc3d2(out_nlpdf.x/phs_scale, N/2, 106, [90,90,-90], [-0.10,0.14], [], 'QSM: wtv-nlpdf')
imagesc3d2(out_nlpdfp.x/phs_scale, N/2, 107, [90,90,-90], [-0.10,0.14], [], 'QSM: wtv-nlpdfp')


imagesc3d2(out_lbv.x/phs_scale-chi_cosmos, N/2, 102, [90,90,-90], [-0.10,0.14], [], 'QSM: wtv-lbv')
imagesc3d2(out_pdf.x/phs_scale-chi_cosmos, N/2, 103, [90,90,-90], [-0.10,0.14], [], 'QSM: wtv-pdf')
imagesc3d2(out_wpdf.x/phs_scale-chi_cosmos, N/2, 104, [90,90,-90], [-0.10,0.14], [], 'QSM: wtv-wpdf')
imagesc3d2(out_wpdfi.x/phs_scale-chi_cosmos, N/2, 105, [90,90,-90], [-0.10,0.14], [], 'QSM: wtv-wpdfi')
%imagesc3d2(out_nlpdf.x/phs_scale, N/2, 106, [90,90,-90], [-0.10,0.14], [], 'QSM: wtv-nlpdf')
imagesc3d2(out_nlpdfp.x/phs_scale-chi_cosmos, N/2, 107, [90,90,-90], [-0.10,0.14], [], 'QSM: wtv-nlpdfp')


% Nonlinear reconstructions
params.alpha1 = 1.0*10^(-(28+20)/10); 
params.mu1 = 1.0*10^(-28/10); 

outnl = nlTV(params);
chi_nl = outnl.x/phs_scale;
metrics_nltv = compute_metrics(real(chi_nl.*mask_use),chi_cosmos);
imagesc3d2(chi_nl, N/2, 4, [90,90,-90], [-0.10,0.14], [], 'QSM: nltv')

outnl2 = nlTGV(params);
chi_nl2 = outnl2.x/phs_scale;
metrics_nltgv = compute_metrics(real(chi_nl2.*mask_use),chi_cosmos);
imagesc3d2(chi_nl2, N/2, 5, [90,90,-90], [-0.10,0.14], [], 'QSM: nltgv')



%% Optimize for RMSE Example


for ka = 1:50
    
    
    alpha = 10^(-ka/10-1);
    params.alpha1 = alpha;
    params.mu1 = 25*alpha;                  % gradient consistency
    
    
params.weight = ones(size(phase_use));
    out = wTV(params); 
    chi = out.x/phs_scale;
    rmse(ka) = compute_rmse(real(chi.*mask_use),chi_cosmos);
    
params.weight = magn_use;
    
    out = wTV(params); 
    chi = out.x/phs_scale;
    rmse_w(ka) = compute_rmse(real(chi.*mask_use),chi_cosmos);
    
    out = nlTV(params); 
    chi = out.x/phs_scale;
    rmse_nl(ka) = compute_rmse(real(chi.*mask_use),chi_cosmos);
end

plot(log(rmse),'b');
hold on;
plot(log(rmse_w),'r');
plot(log(rmse_nl),'k');
hold off;

ka = 5;
params.weight = ones(size(phase_use));
    alpha = 10^(-ka/10-1);
    params.alpha1 = alpha;
    params.mu1 = 100*alpha;                  % gradient consistency
    out = wTV(params); 
    chi = out.x/phs_scale;
    
ka = 27;
params.weight = magn_use;
    alpha = 10^(-ka/10-1);
    params.alpha1 = alpha;
    params.mu1 = 100*alpha;                  % gradient consistency
    outw = wTV(params); 
    chiw = outw.x/phs_scale;
    outnl = wTV(params); 
    chinl = outnl.x/phs_scale;
    
    
    
    
    %%
    
params.beta = 10.^( ((19)-15)/4 +2);
for ss = 1:30
    
    
%params.beta = 10.^( ((ss)-15)/4 +2);
params.muh = 10.^( ((ss)-15)/4 +2);%1000*params.beta;
outss = sswTV(params); 
    chi = outss.x/phs_scale;
    
    rmse_ss(ss) = compute_rmse(real(chi.*mask_use),chi_cosmos);
end

params.beta = 10.^( ((19)-15)/4 +2); %1000
params.muh = 10.^( ((26)-15)/4 +2); %56234. Ratio = 56.2

outss = sswTV(params); 
    chi = outss.x/phs_scale;
    metrics_ss = compute_metrics(real(chi.*mask_use),chi_cosmos);
imagesc3d2(chi, N/2, 6, [90,90,-90], [-0.10,0.14], [], 'QSM: sstv')
imagesc3d2(outss.phi, N/2, 7, [90,90,-90], [-0.10,0.14], [], 'QSM: phi')
imagesc3d2(chi-out.x/phs_scale, N/2, 8, [90,90,-90], [-0.10,0.14], [], 'QSM: diff')