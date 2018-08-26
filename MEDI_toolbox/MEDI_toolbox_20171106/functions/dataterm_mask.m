% Generate the data weighting
%   w=dataterm_mask(dataterm_weighting_mode, N_std, Mask,cutoff)
% 
%   output
%   w - data weighting 
%
%   input
%   dataterm_weighting_mode - 0, uniform weighting; 1, SNR weighting
%   N_std - noise standard deviation on the field map
%   Mask is a binary 3D matrix denoting the Region Of Interest
%
%   Created by Ildar Khalidov in 20010
%   Last modified by Tian Liu on 2013.07.24


function w=dataterm_mask(dataterm_weighting_mode, N_std, Mask)


switch dataterm_weighting_mode
case 0
    w=1;
case 1 
    w = Mask./N_std;
    w(isnan(w)) = 0;
    w(isinf(w)) = 0;
    w = w.*(Mask>0);
    w = w./mean(w(Mask>0)); %normalize w by its mean    

end
