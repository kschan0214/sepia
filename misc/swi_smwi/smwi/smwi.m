%% [pSMWI, dSMWI] = smwi(magn,qsm,thres,m)
%
% Input
% --------------
% magn          : image that combines with QSM contrast, can be 3D/4D
% qsm           : QSM map (in ppm)
% thres         : threshold used to create para-/dia-magnetic contrast
% m             : power to enhance QSM contrast
%
% Output
% --------------
% pSMWI         : paramagnetic susceptibility map weighted image
% dSMWI         : diamagnetic susceptibility map weighted image
%
% Description: compute susceptibility map-weighted imaging
% Assume data dimension: [x y z t]
% if no input of thres and m, recommended thres and m will be used
% ref: MRM 72:337-346(2014)
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 7 October 2016
% Date modified: 14 April 2019
%
%
function [pSMWI, dSMWI] = smwi(magn,qsm,thres,m)

magn = double(magn);
qsm  = double(qsm);

if nargin < 4
    m = 4; %power
end
if nargin < 3 || thres < 0
    thres = 1; %ppm
    warning('Threshold has to be greater than zero. Using 1 ppm instead...');
end

% display algorithm parameters
disp(['Threshold (ppm): ' num2str(thres)]);
disp(['Contrast power: ' num2str(m)]);

% parammagnetic mask
pSmask              = (thres-qsm)/thres;
pSmask(qsm>thres)   = 0; 
pSmask(qsm<0)       = 1;

% diamagnetic mask
dSmask              = (thres+qsm)/thres;
dSmask(qsm>0)       = 1; 
dSmask(qsm<-thres)  = 0;


pSMWI = bsxfun(@times,magn,pSmask.^m);
dSMWI = bsxfun(@times,magn,dSmask.^m);

end
