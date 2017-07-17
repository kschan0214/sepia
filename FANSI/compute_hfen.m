% Calculates the High Frequency Error metric.
%
% See README.txt for further information.
% Last modified by Carlos Milovic in 2017.03.30
%
function [ hfen ] = compute_hfen( img1, img2 )
 

% Laplacian of Gaussian filter to get high frequency information:

filt_siz = [1,1,1] * 15;
sig = [1,1,1] * 1.5;

siz = (filt_siz-1)/2;
[x,y,z] = ndgrid(-siz(1):siz(1), -siz(2):siz(2), -siz(3):siz(3));

h = exp(-(x.*x/2/sig(1)^2 + y.*y/2/sig(2)^2 + z.*z/2/sig(3)^2));
h = h / sum(h(:));

arg = (x.*x/sig(1)^4 + y.*y/sig(2)^4 + z.*z/sig(3)^4 - (1/sig(1)^2 + 1/sig(2)^2 + 1/sig(3)^2));
H = arg .* h;
H = H - sum(H(:)) / prod(2*siz+1);
        
        
img1_log = imfilter(img1, H, 'same');

img2_log = imfilter(img2, H, 'same');


hfen = compute_rmse(img1_log, img2_log);


end

