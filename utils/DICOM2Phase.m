%% phase = DICOM2Phase(dicomPhase)
%
% Usage:
%
% Input
% --------------
%
% Output
% --------------
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 11 April 2018
% Date last modified:13 April 2-18
%
%
function phase = DICOM2Phase(niiPhase)
dicomPhase = double(niiPhase.img);
scaleSlope=niiPhase.hdr.dime.scl_slope;
scaleIntercept=niiPhase.hdr.dime.scl_inter;

newMax = niiPhase.hdr.dime.glmax*2 + scaleIntercept;
newMin = niiPhase.hdr.dime.glmin*2 + scaleIntercept;
fullRange = newMax-newMin + 1;

% scale to true value of nifti file
dicomPhaseRescale = (dicomPhase*scaleSlope) + scaleIntercept ;
% Scale the full range to [-pi,pi)
phase = (dicomPhaseRescale-scaleIntercept) / fullRange * 2*pi - pi;

end