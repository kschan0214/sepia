%% BipolarCorr = BipolarEddyCorrect(Bipolar,unwrap)
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
% Original: EddyCurrentVisualisationCorrection.m from Jose P. Marques
% Modified by Kwok-Shing Chan
% k.chan@donders.ru.nl
% Date created: 13 April 2018
% Date last modified: 5 June 2019
%
%
function bipolarCorr = BipolarEddyCorrect(bipolarCplxME,mask,algorParam)

sepia_universal_variables;

algorParam = check_and_set_SEPIA_algorithm_default(algorParam);

disp('--------------------------');
disp('Bipolar readout correction');
disp('--------------------------');
disp('Correcting eddy current effect on bipolar readout data...');

bipolarCplxME(isnan(bipolarCplxME))=0;

dims = size(bipolarCplxME) ;

% successive phase difference mean made so that the last echo taken into
% account has to be an odd number ()
% to= dims(4)-(-mod(dims(4),2)+1);
% Phase of mean evens minus mean odds - contains Eddy Current and
to = dims(4)-(mod(dims(4),2));
PhaseDiffEvensMinusOdds = angle(mean(bipolarCplxME(:,:,:,2:2:to)./bipolarCplxME(:,:,:,1:2:to-1),4));
PhaseDiffEvensMinusOdds(isnan(PhaseDiffEvensMinusOdds)) = 0;

% Phase of even minus odd images - contains Eddy Current and Field Map
% PhaseDiffEvensMinusOddsFluctuations=angle(mean(bipolarCplxME(:,:,:,2:2:end-2)./bipolarCplxME(:,:,:,1:2:end-3),4));

PhaseDiffEvens  = angle(mean(bipolarCplxME(:,:,:,4:2:end)./bipolarCplxME(:,:,:,2:2:end-2),4));
PhaseDiffOdds   = angle(mean(bipolarCplxME(:,:,:,3:2:end)./bipolarCplxME(:,:,:,1:2:end-2),4));
PhaseDiffEvens(isnan(PhaseDiffEvens))=0;
PhaseDiffOdds(isnan(PhaseDiffOdds))=0;
% BipolarOddvsEven=mean(abs(bipolarCplxME(:,:,:,1:2:end)),4)-mean(abs(bipolarCplxME(:,:,:,2:2:end-1)),4);
% this has no unwrapping done to it - but it is a good first estimate
% meanEddycurrent = PhaseDiffEvensMinusOdds-0.25*(PhaseDiffEvens+PhaseDiffOdds);

%% a more robust eddy current estimation has to do unwrapping
te=1:dims(4);

algorParam.unwrap.unit 	= 'radHz';
algorParam.unwrap.echoCombMethod = methodEchoCombineName{1}; 

headerAndExtraData.te   = te(4)-te(2);
headerAndExtraData.magn	= double(mean(abs(bipolarCplxME(:,:,:,2:2:end-2)),4));
[fieldMapEven,~] = estimateTotalField(double(PhaseDiffEvens),mask,dims(1:3),[1 1 1],algorParam,headerAndExtraData);

headerAndExtraData.te   = te(3)-te(1);
headerAndExtraData.magn	= double(mean(abs(bipolarCplxME(:,:,:,1:2:end-2)),4));
[fieldMapOdd,~] = estimateTotalField(double(PhaseDiffOdds),mask,dims(1:3),[1 1 1],algorParam,headerAndExtraData);

% [fieldMapEven,~] = estimateTotalField(double(PhaseDiffEvens),double(mean(abs(bipolarCplxME(:,:,:,2:2:end-2)),4)),dims(1:3),[1 1 1],...
%                         'Unwrap',unwrap,'TE',te(4)-te(2),'unit','radHz','mask',mask);
% [fieldMapOdd,~] = estimateTotalField(double(PhaseDiffOdds),double(mean(abs(bipolarCplxME(:,:,:,1:2:end-2)),4)),dims(1:3),[1 1 1],...
%                         'Unwrap',unwrap,'TE',te(3)-te(1),'unit','radHz','mask',mask);


% Even=cat(4,mean(abs(bipolarCplxME(:,:,:,2:2:end-2)),4),mean(abs(bipolarCplxME(:,:,:,4:2:end)),4).*exp(i*PhaseDiffEvens)) ;
% Odd=cat(4,mean(abs(bipolarCplxME(:,:,:,1:2:end-2)),4),mean(abs(bipolarCplxME(:,:,:,3:2:end)),4).*exp(i*PhaseDiffOdds)) ;
% 
% [fieldMapEven,~] = estimateTotalField(double(angle(Even)),double(abs(Even)),dims(1:3),[1 1 1],...
%                         'Unwrap',unwrap,'TE',te(2:2:4),'unit','radHz','mask',mask);
% [fieldMapOdd,~] = estimateTotalField(double(angle(Odd)),double(abs(Odd)),dims(1:3),[1 1 1],...
%                         'Unwrap',unwrap,'TE',te(1:2:3),'unit','radHz','mask',mask);

% [fieldMapEven, ~, ~, ConfidenceEven]=T2starAndFieldCalc(abs(Even),angle(Even),te(2:2:4));
% [fieldMapOdd, ~, ~, ConfidenceOdd]=T2starAndFieldCalc(abs(Odd),angle(Odd),te(1:2:3));
% % [fieldMapEven, ~, ~, ConfidenceEven]=T2starAndFieldCalc(abs(Bipolar(:,:,:,2:2:end)),angle(Bipolar(:,:,:,2:2:end)),te(2:2:end));
% % [fieldMapOdd, ~, ~, ConfidenceOdd]=T2starAndFieldCalc(abs(Bipolar(:,:,:,1:2:end)),angle(Bipolar(:,:,:,1:2:end)),te(1:2:end));


% the last echo to be taken into account has to be even
headerAndExtraData.te   = te(1:2);
headerAndExtraData.magn	= double(mean(abs(bipolarCplxME(:,:,:,1:2:to)),4));
[fieldMapOddEven,~] = estimateTotalField(double(PhaseDiffEvensMinusOdds),mask,dims(1:3),[1 1 1],algorParam,headerAndExtraData);

% [fieldMapOddEven,~] = estimateTotalField(double(PhaseDiffEvensMinusOdds),double(mean(abs(bipolarCplxME(:,:,:,1:2:to)),4)),dims(1:3),[1 1 1],...
%                         'Unwrap',unwrap,'TE',te(1:2),'unit','radHz','mask',mask);

% OddEven=cat(4,mean(abs(bipolarCplxME(:,:,:,1:2:to)),4),mean(abs(bipolarCplxME(:,:,:,2:2:to)),4).*exp(i*PhaseDiffEvensMinusOdds)) ;
% 
% [fieldMapOddEven,~] = estimateTotalField(double(angle(OddEven)),double(abs(OddEven)),dims(1:3),[1 1 1],...
%                         'Unwrap',unwrap,'TE',te(1:2),'unit','radHz','mask',mask);
% [fieldMapOddEven, ~, ~, ConfidenceOdd]=T2starAndFieldCalc(abs(OddEven),angle(OddEven),te(1:2));

meanEddycurrent2 = fieldMapOddEven-0.5*(fieldMapOdd+fieldMapEven);

% EddycurrentVariations=angle(...
%     bsxfun(@times,bsxfun(@times,bipolarCplxME(:,:,:,2:2:end-1) , ...
%     2./(bipolarCplxME(:,:,:,1:2:(end-2))+bipolarCplxME(:,:,:,3:2:end))),exp(-i * meanEddycurrent2))...
%     );

%%

% the chosen fitting strategy takes the magnitude image as weights... and
% nothing else fancy.. first order correction gave the best results
meanEddycurrent2(isnan(meanEddycurrent2))=0;

% [FIT3D,~,~]=PolyFit(meanEddycurrent2,abs(sum(bipolarCplxME,4)),1);
% [FIT3D,~,~]=PolyFit(meanEddycurrent2,sum(abs(bipolarCplxME),4),1);
[FIT3D,~,~]=PolyFit(double(meanEddycurrent2),mask,1);


bipolarCorr = zeros(dims);
for echo=1:dims(4)
    bipolarCorr(:,:,:,echo)=bipolarCplxME(:,:,:,echo).*exp(-1i *(-1)^echo * 0.5 * FIT3D);
end

end
