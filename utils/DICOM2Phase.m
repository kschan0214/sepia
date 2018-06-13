%% phase = DICOM2Phase(niiPhase)
%
% Input
% --------------
% niiPhase      : NIfTI matlab structure contains phase image directly
%                 converted from DICOM
%
% Output
% --------------
% phase         : phase image rescaled to [-pi,pi)
%
% Description: This function convert the DICOM phase values (usually
% [-4096,4095]) to radian ([-pi,pi)])
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 11 April 2018
% Date last modified: 12 June 2018
%
%
function phase = DICOM2Phase(niiPhase)
% change datatype to double
dicomPhase = double(niiPhase.img);
% get scale slope
scaleSlope=niiPhase.hdr.dime.scl_slope;
% get scale inetercept
scaleIntercept=niiPhase.hdr.dime.scl_inter;

% compute the new maximum value
newMax = niiPhase.hdr.dime.glmax*scaleSlope + scaleIntercept;
% compute the new minimum value
newMin = niiPhase.hdr.dime.glmin*scaleSlope + scaleIntercept;
% +1 for 0
fullRange = newMax-newMin + 1;

% scale to true value of nifti file
dicomPhaseRescale = (dicomPhase*scaleSlope) + scaleIntercept ;
% scale the full range to [-pi,pi)
phase = (dicomPhaseRescale-newMin) / fullRange * 2*pi - pi;

end