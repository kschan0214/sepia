function [ y ] = SSQSM_vsharp( X, D, Del_Sharp, cDel_Sharp, DiffMask, mu, fdx, fdy, fdz, cfdx, cfdy, cfdz, magn_weight )
%APPLY_FORWARD Summary of this function goes here
%   Detailed explanation goes here

x = reshape(X, size(D));
Dx = D .* x;

tmp = 0;
for k = 1:size(Del_Sharp,4)
    tmp = tmp + DiffMask(:,:,:,k) .* ifftn(Del_Sharp(:,:,:,k) .* Dx);
end


temp = 0;
for k = 1:size(Del_Sharp,4)
    temp = temp + cDel_Sharp(:,:,:,k) .* fftn(DiffMask(:,:,:,k) .* tmp);
end
temp = conj(D) .* temp;


y = temp + mu .* (cfdx.*fftn(magn_weight(:,:,:,1) .* ifftn(fdx.*x)) + ...
    cfdy.*fftn(magn_weight(:,:,:,2) .* ifftn(fdy.*x)) + cfdz.*fftn(magn_weight(:,:,:,3) .* ifftn(fdz.*x)));

y = y(:);


end
