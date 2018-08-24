% Noise standard deviation calculation
%   noise_level = calfieldnoise(iField, Mask)
% 
%   output
%   noise_level - the noise standard devivation of complex MR signal
%
%   input
%   iField - the complex MR image
%   Mask - the a region of interest that is not pure noise
%
%   Created by Tian Liu in 20013.07.24
%   Last modified by Tian Liu on 2013.07.24

function noise_level = calfieldnoise(iField, Mask)
expected_SNR = 40;
iMag = sqrt(sum(abs(iField).^2,4));
iField1 = iField(:,:,:,1).*(~Mask).*(iMag<max(iMag(:))/expected_SNR).*(iMag>0);
noise_level = std(real(iField1(iField1~=0))); 
