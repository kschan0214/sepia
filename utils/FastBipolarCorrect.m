%% BipolarCorr = FastBipolarCorrect(Bipolar,unwrap)
%
% Input
% --------------
% bipolarCplxME : complex-valued multi-echo GRE data with bipolar readout
% unwrap        : phase unwrapped method (optional)
%
% Output
% --------------
% bipolarCorr   : eddy current corrected data
%
% Description: correct the inconsistent phase between odd and even echoes
% due to bipolar greadout
%
% Original: FastBipolarCorrect.m from Jose P. Marques
% Modified by Kwok-Shing Chan
% k.chan@donders.ru.nl
% Date created: 5 August 2022
% Date modified: 
%
%
function [img_cplx_corr,FIT3D_weights]  = FastBipolarCorrect(img_cplx,mask)

% get image size
dims = size(img_cplx);
% This is to ensure the same number of even and odd echos are used in the estimation
to = dims(4)-(mod(dims(4),2));
% compute phase difference between odd and even echoes
% PhaseDiffEvensMinusOdds has two contributions: (1) bipolar readout
% induced phase (phi_bi) and (2) susceptibility induced phase  (phi_x), in other
% words, PhaseDiffEvensMinusOdds = phi_bi + phi_x
PhaseDiffEvensMinusOdds = angle(mean(img_cplx(:,:,:,2:2:to)./img_cplx(:,:,:,1:2:to-1),4));
PhaseDiffEvensMinusOdds(isnan(PhaseDiffEvensMinusOdds)) = 0;

% compute the phase induced by susceptibility, by first separating odd and
% even echoes and then combine the results to boost SNR
% PhaseDiff = 2*phi_x
PhaseDiff  = angle(sum(img_cplx(:,:,:,4:2:end)./img_cplx(:,:,:,2:2:end-2),4) + sum(img_cplx(:,:,:,3:2:end)./img_cplx(:,:,:,1:2:end-2),4));
PhaseDiff(isnan(PhaseDiff)) = 0;

% compute bipolar readout induced phase
% meanEddycurrent2 = 2*PhaseDiffEvensMinusOdds - PhaseDiff = 2*phi_bi + 2*phi_x - 2*phi_x = 2*phi_bi
% meanEddycurrent = fieldMapOddEven-0.5*(fieldMapOdd+fieldMapEven);
% TODO 20220912: phase wrapping can happen in meanEddycurrent2
meanEddycurrent2 = angle(exp(1i * (PhaseDiffEvensMinusOdds + PhaseDiffEvensMinusOdds - PhaseDiff)));
meanEddycurrent2(isnan(meanEddycurrent2))=0;

% uses a weigthing based on the mean signal form the 3rd echo onwards
%
weights = mean(abs(img_cplx(:,:,:,3:end)),4);
weights (weights>prctile(weights(:),90))=prctile(weights(:),90);
% [FIT3D_weights,~,b2]=PolyFit_weights(double(meanEddycurrent2)/2,mask,1,double(weights));

% estimate the 1st order coefficients first
fit3D_1st_init = first_order_polynomial_coef(meanEddycurrent2,mask);
% remove the 1st order coreffients from meanEddycurrent2 to avoid phase jumps in Polyfit
meanEddycurrent2_0th = angle(exp(1i*meanEddycurrent2) ./ exp(1i*fit3D_1st_init));
% remove the mean offset to avoid phase jump 
tmp = exp(1i*meanEddycurrent2_0th) ./ exp(1i*mean(meanEddycurrent2_0th( mask == 1)));
% further refine the results by perfoming 1st order PolyFit, now the
% residual tmp should have less phase wrapping issue
[FIT3D,~,~] = PolyFit_weights(angle(tmp),mask,1,double(weights));
% [FIT3D,~,~] = PolyFit(angle(tmp2),mask,1);
FIT3D_weights = (fit3D_1st_init + mean(meanEddycurrent2_0th( mask == 1)) + FIT3D)/2;
% FIT3D_weights = (fit3D_1st_init + FIT3D)/2;

compare_to_other_methods = 0;
if compare_to_other_methods==1
    [FIT3D_old,~,b1]=PolyFit_weights(meanEddycurrent2/2,mask,1,[]);
    [FIT3D_maskFree,~,b3]=PolyFit_weights(double(meanEddycurrent2/2),ones(dims(1:3)),1,double(weights));

    subplot(322)
    Orthoview(abs(FIT3D_old-meanEddycurrent2/2),[],[0 0.1]/2);
    title('Residual linear fit inside mask')
    subplot(324)
    Orthoview(abs(FIT3D_weights-meanEddycurrent2/2),[],[0 0.1]/2);
    title(' Residual of linear fit inside mask using weighting function')
    subplot(323)
    Orthoview(weights.*mask,[],[]);
    title('weight times mask')
    subplot(326)
    Orthoview(abs(FIT3D_maskFree-meanEddycurrent2/2),[],[0 0.1]/2);
    title('linear fit all volume')
    subplot(325)
    Orthoview(weights,[],[]);
    title('weightwithout mask')
end;

img_cplx_corr = zeros(size(img_cplx));
for echo=1:dims(4)
    img_cplx_corr(:,:,:,echo)=img_cplx(:,:,:,echo).*exp(-1i *(-1)^echo * 0.5 * (FIT3D_weights));
end

PhaseDiffEvensMinusOdds = mean(img_cplx_corr(:,:,:,2:2:to)./img_cplx_corr(:,:,:,1:2:to-1),4);

PhaseDiffEvensMinusOdds_mean = mean(PhaseDiffEvensMinusOdds(and(and(mask==1,~isnan(PhaseDiffEvensMinusOdds)),~isinf(PhaseDiffEvensMinusOdds))));

if abs(angle(PhaseDiffEvensMinusOdds_mean ./ exp(1i*pi))) < abs(angle(PhaseDiffEvensMinusOdds_mean))
    FIT3D_weights = FIT3D_weights - pi;
    for echo=1:dims(4)
        img_cplx_corr(:,:,:,echo)=img_cplx(:,:,:,echo).*exp(-1i *(-1)^echo * 0.5 * (FIT3D_weights));
    end
end
    
FIT3D_weights = FIT3D_weights .* mask;

end

function fit3D_1st = first_order_polynomial_coef(phase_wrapped,mask)

b = zeros(4,1);

% convert phase into complex valued
phase_wrapped = exp(1i*phase_wrapped);

% compute coefficients along three directions
for k = 1:3
    % compute phase difference (gradient, 1st order Polynomial) along a single dimension
    phase_diff = angle(phase_wrapped ./ (circshift(phase_wrapped,1,k)));
    
    % fit a constant (0-th order) on gradient -> 1st order coefficient
%     b(k) = PolyFit_weights(double(phase_diff),mask,0,double(weights));
%     [~,~,b(k+1)] = PolyFit(double(phase_diff),mask,0);
    b(k+1) = mean(phase_diff(mask==1));
%     [~,~,b(k+1)] = PolyFit_weights(double(phase_diff),mask,0,[]);
end
dim             = size(mask);
Indices         = find(ones(dim));
[x1,y1,z1]      = ind2sub(dim,Indices);
model_1storder  = Create_1st_order_model(x1,y1,z1,dim);
Fit             = model_1storder*b;
fit3D_1st    	= reshape(Fit,dim);

end

function model = Create_1st_order_model(x1,y1,z1,dim)
model = double(zeros(length(x1),4));
% 0th order coeffient
model(1:length(x1),1)=1;
% 1st order coeffients
model(:,2)=reshape(x1-dim(1)/2,length(x1),1);%x -siemens
model(:,3)=reshape(y1-dim(2)/2,length(x1),1);%y -siemens
model(:,4)=reshape(z1-dim(3)/2,length(x1),1);%z -siemens
end