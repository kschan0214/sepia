function [img_high, img_low, FW] = CalcPhaseHanning(img,hanning_Diameter)
% Wei Li, Duke University, 6/5/2010
img=padarray(img,[0 0 71]);
SS = size(img);
yy=[-SS(2)/2:SS(2)/2-1]/SS(2);
xx=[-SS(1)/2:SS(1)/2-1]/SS(1);
zz=[-SS(3)/2:SS(3)/2-1]/SS(3);
[y,x,z] = meshgrid(yy,xx,zz);
krad = sqrt(x.^2+y.^2+z.^2);
N=hanning_Diameter;
FW=zeros(SS);
AA=0.5*(1+cos(2*pi*krad/N));
FW(krad<(N/2))=AA(krad<(N/2));
img_low = ifftnc(fftnc(img(:,:,:)).*FW);
img_high=img-img_low;
img_low=img_low(:,:,72:end-71);
img_high=img_high(:,:,72:end-71);

