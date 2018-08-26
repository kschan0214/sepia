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

load spatial_res;           % voxel size
load msk;                   % brain mask => obtained by eroding the BET mask by 5 voxels (by setting peel=5 in LBV)


load magn;                  % magnitude from transversal orientation
load chi_cosmos;            % COSMOS from 12 orientations (in ppm)   
load ext_masks;            % COSMOS from 12 orientations (in ppm)   

N = size(chi_cosmos);

center = N/2 + 1;

TE = 25e-3;
B0 = 2.8936;
gyro = 2*pi*42.58;

phs_scale = TE * gyro * B0;

imagesc3d2(chi_cosmos, N/2, 1, [90,90,-90], [-0.10,0.14], [], 'COSMOS')

chi_full = chi_cosmos;
chi_full((1-mask_use)>0) = 9.395/1;
chi_full(msk_tissue>0) = 0;
chi_full(msk_deleted>0) = 0;
imagesc3d2(chi_full, N/2, 1, [90,90,-90], [-0.50,1.0], [], 'COSMOS')


%% Create dipole kernel and susceptibility to field model

kernel = dipole_kernel( N, spatial_res, 0 ); % 0 for the continuous kernel by Salomir and Marques and Bowtell.
                                             % 1 for a discrete kernel
                                             % formulation
                                             % 2 for an integrated Green
                                             % function


% Create raw phase
rphase = phs_scale*ifftn(fftn(chi_full).*kernel);
mp = mean(rphase(:));
lphase = phs_scale*ifftn(fftn(chi_cosmos).*kernel);
bphase = phs_scale*ifftn(fftn(chi_full-chi_cosmos).*kernel);

signal = magn.*exp(1i*(rphase-mp))+0.025*(randn(size(magn))+1i*randn(size(magn))); %SNR = 40
phase_use = angle(signal);
magn_use = magn/max(magn(:));
mask_use = msk;

imagesc3d2(phase_use, N/2, 10, [90,90,-90], [-3.14,3.14], [], 'Phase')
imagesc3d2(angle(Img), N/2, 11, [90,90,-90], [-3,3], [], 'Phase single acq')
imagesc3d2(mask_use, N/2, 3, [90,90,-90], [0,1], [], 'Mask')
imagesc3d2(mask_use+msk_deleted, N/2, 4, [90,90,-90], [0,1], [], 'Mask')
imagesc3d2(magn_use, N/2, 3, [90,90,-90], [0,1], [], 'Magn')

%% Unwrapping results
phase_unwrap = unwrap(phase_use,spatial_res);
phase_unwrapL = unwrapLaplacian(phase_use,N,spatial_res);


imagesc3d2((rphase-mp).*mask_use, N/2, 20, [90,90,-90], [-6,6], [], 'Phase')
imagesc3d2(phase_unwrapL.*mask_use, N/2, 21, [90,90,-90], [-6,6], [], 'PuwL')
imagesc3d2(phase_unwrap.*mask_use, N/2, 22, [90,90,-90], [-6,6], [], 'Puw')

phase_Munwrap = unwrap(phase_use.*mask_use,spatial_res);
phase_MunwrapL = unwrapLaplacian(phase_use.*mask_use,N,spatial_res);

imagesc3d2(phase_MunwrapL.*mask_use, N/2, 23, [90,90,-90], [-6,6], [], 'PuwL')
imagesc3d2(phase_Munwrap.*mask_use, N/2, 24, [90,90,-90], [-6,6], [], 'Puw')

% Unwrapped filtering
nlocalp = (phase_unwrap-bphase-mean(phase_unwrap(:))+mean(bphase(:))).*mask_use;
nlocalp(nlocalp>pi) = nlocalp(nlocalp>pi)-2*pi;
nlocalp(nlocalp<-pi) = nlocalp(nlocalp<-pi)+2*pi;
nlocalp = (nlocalp-sum(nlocalp(:))/sum(mask_use(:))+sum(lphase(:))/sum(mask_use(:))).*mask_use;
imagesc3d2(nlocalp, N/2, 29, [90,90,-90], [-1.12,1.12], 0, 'local phase');
imagesc3d2(lphase, N/2, 30, [90,90,-90], [-1.12,1.12], 0, 'local phase');

pLBV = LBV(phase_unwrapL,mask_use,N,spatial_res);
imagesc3d2(pLBV, N/2, 31, [90,90,-90], [-1.12,1.12], 0, 'pLBV');
imagesc3d2(pLBV-nlocalp, N/2, 31, [90,90,-90], [-1.12,1.12], 0, 'pLBV');

[lpPDF hpPDF] = PDF(phase_unwrapL, 1./(40*magn_use+eps), mask_use,N,spatial_res, [0 0 1]);
imagesc3d2(lpPDF, N/2, 32, [90,90,-90], [-1.12,1.12], 0, 'PDF');
imagesc3d2(lpPDF-nlocalp, N/2, 32, [90,90,-90], [-1.12,1.12], 0, 'PDF');


pLBVi = LBV(phase_Munwrap,mask_use,N,spatial_res);
imagesc3d2(pLBVi, N/2, 33, [90,90,-90], [-1.12,1.12], 0, 'pLBV');
imagesc3d2(pLBVi-nlocalp, N/2, 33, [90,90,-90], [-1.12,1.12], 0, 'pLBV');

[lpPDFi hpPDFi] = PDF(phase_Munwrap, 1./(40*magn_use+eps), mask_use,N,spatial_res, [0 0 1]);
imagesc3d2(lpPDFi, N/2, 34, [90,90,-90], [-1.12,1.12], 0, 'PDF');
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
ppdf.maxOuterIter = 150;
ppdf.tol_update = 0.01;
ppdf.K = kernel;
ppdf.input = phase_use;
ppdf.mask = mask_use;
ppdf.outermask = mask_use+msk_deleted+msk_tissue;
ppdf.weight = magn_use.*mask_use;
ppdf.mu1 = 0.0018;                  % gradient consistency
ppdf.alpha1 = eps;%0.0001;              % gradient L1 penalty
nout5 = nlPDFp(ppdf);

imagesc3d2(nout5.local.*mask_use, N/2, 71, [90,90,-90], [-1.12,1.12], 0, 'nlPDF');
imagesc3d2(nout5.local.*mask_use-nlocalp, N/2, 71, [90,90,-90], [-1.12,1.12], 0, 'nlPDF');
nout5t = nlPDFp2(ppdf);

imagesc3d2(nout5.phi.*mask_use, N/2, 72, [90,90,-90], [-6,6], 0, 'nlPDF');

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