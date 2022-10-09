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
% Description: This function converts the DICOM phase values (usually
% [-4096,4095]) to radian ([-pi,pi)]), assuming the maximum and minimum
% value in the data corresponding to pi and pi.
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 11 April 2018
% Date modified: 12 June 2018
% Date modified: 4 August 2022 (v1.0.1)
%
%
function phase = DICOM2Phase(niiPhase)
% % change datatype to double
% dicomPhase = double(niiPhase.img);
% get scale slope
scaleSlope=niiPhase.hdr.dime.scl_slope;
% get scale inetercept
scaleIntercept=niiPhase.hdr.dime.scl_inter;

% rescale phase
phase_rescale = single(niiPhase.img) * scaleSlope + scaleIntercept;
% find maximum and minimum in the data after rescaling
max_rescale = max(phase_rescale(:));
min_rescale = min(phase_rescale(:));

% rescale phase data between -pi and pi usng the max and min in the data
% the limitation of this method is that if the max and min values might not
% corresponding to pi and -pi
phase = (phase_rescale - min_rescale) / (max_rescale - min_rescale) *2*pi - pi;

% % compute the new maximum value
% newMax = niiPhase.hdr.dime.glmax*scaleSlope + scaleIntercept;
% % compute the new minimum value
% newMin = niiPhase.hdr.dime.glmin*scaleSlope + scaleIntercept;
% % +1 for 0
% fullRange = newMax-newMin + 1;
% 
% % scale to true value of nifti file
% dicomPhaseRescale = (dicomPhase*scaleSlope) + scaleIntercept ;
% % scale the full range to [-pi,pi)
% phase = (dicomPhaseRescale-newMin) / fullRange * 2*pi - pi;

end