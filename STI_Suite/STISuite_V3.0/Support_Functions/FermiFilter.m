function [img_high, img_low,FW] = FermiFilter(img,filterRatio)
% Fermi Filtering
% Wei Li, Duke University, 6/5/2010

tic
SS = size(img);
[y,x,z] = meshgrid(-SS(2)/2:SS(2)/2-1,-SS(1)/2:SS(1)/2-1,-SS(3)/2:SS(3)/2-1);
fermir=filterRatio*SS(1)/2;
fermiw=0.1*fermir;
krad = sqrt(x.^2+y.^2+z.^2);
FW = 1./(1+exp((krad-fermir)/fermiw));
img_low = ifftnc(fftnc(img(:,:,:)).*FW);
img_high= img-img_low;
toc
