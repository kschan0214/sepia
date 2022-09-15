function ss  = FastFirstEchoEddyCorrect(ss,mask);

dims = size(ss);
to = dims(4)-(mod(dims(4),2)); % This is to ensure the same number of even and odd echos are used in the estimation
PhaseDiff_2_1 = (abs(ss(:,:,:,2)).*(ss(:,:,:,2)./ss(:,:,:,1)));
PhaseDiff_2_1(isnan(PhaseDiff_2_1)) = 0;

PhaseDiff_end_2  = (mean(abs(ss(:,:,:,3:end)).*exp(i*angle(ss(:,:,:,3:end)./ss(:,:,:,2:end-1))),4));
PhaseDiff_end_2(isnan(PhaseDiff_end_2)) = 0;


% phase_diff_2end_vs_1 = angle(mean(abs(ss(:,:,:,3:end)).*exp(i*angle(ss(:,:,:,3:end)./ss(:,:,:,2:end-1))),4)./...
%     (ss(:,:,:,2)./ss(:,:,:,1)));
% 

meanEddycurrent2 = angle(PhaseDiff_end_2./PhaseDiff_2_1);
meanEddycurrent2(isnan(meanEddycurrent2))=0;
% uses a weigthing based on the mean signal form the 3rd echo onwards
%
weights = mean(abs(ss(:,:,:,3:end)),4);
weights (weights>prctile(weights(:),90))=prctile(weights(:),90);
[FIT3D_weights,~,b2]=PolyFit_weights(double(meanEddycurrent2),mask,1,double(weights));

compare_to_other_methods = 0;
if compare_to_other_methods==1
    [FIT3D_old,~,b1]=PolyFit_weights(meanEddycurrent2,mask,1,[]);
    [FIT3D_maskFree,~,b3]=PolyFit_weights(double(meanEddycurrent2),ones(dims(1:3)),1,double(weights));

    subplot(322)
    Orthoview(abs(FIT3D_old-meanEddycurrent2),[],[0 0.1]/2);
    title('Residual linear fit inside mask')
    subplot(324)
    Orthoview(abs(FIT3D_weights-meanEddycurrent2),[],[0 0.1]/2);
    title(' Residual of linear fit inside mask using weighting function')
    subplot(323)
    Orthoview(weights.*mask,[],[]);
    title('weight times mask')
    subplot(326)
    Orthoview(abs(FIT3D_maskFree-meanEddycurrent2),[],[0 0.1]/2);
    title('linear fit all volume')
    subplot(325)
    Orthoview(weights,[],[]);
    title('weightwithout mask')
end;

for echo=2:dims(4)
    ss(:,:,:,echo)=ss(:,:,:,echo).*exp(1i * FIT3D_weights);
end
