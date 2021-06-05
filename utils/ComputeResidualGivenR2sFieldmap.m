function residual = ComputeResidualGivenR2sFieldmap(te,r2s,fieldmap,s_cplx)

magn = abs(s_cplx);
phase = angle(s_cplx);

% estimated T1w signal S0
s0 = s_cplx(:,:,:,1).*exp(r2s*te(1));

% compute simulated signal
shat = zeros(size(magn));
for kt=1:length(te)
    shat(:,:,:,kt) = s0.*exp(te(kt)*(-r2s+1i*2*pi*fieldmap));
end

% simulated signal, get rid of the initial phase term
shat = bsxfun(@times,shat,conj(shat(:,:,:,1)));

% measurement, get rid of the initial phase term
s = bsxfun(@times,magn.*exp(1i*phase),conj(magn(:,:,:,1).*exp(1i*phase(:,:,:,1))));

%% compute relative residual
% this formulation emphasises the difference of each echo
residual=sum(abs(shat-s).^2,4)./sum(abs(s).^2,4);  
% this formulation emphasises the Guassian noise of the data
% residual=(abs(sum(shat-s,4)).^2)./sum(abs(s).^2,4);  

residual(isnan(residual)) = max(residual(:));
residual(isinf(residual)) = max(residual(:));

end