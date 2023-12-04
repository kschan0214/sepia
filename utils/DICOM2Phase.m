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
% Date modified: 3 December 2023 (v1.2.6.6)
%
function phase = DICOM2Phase(niiPhase)

% suggested by Simon
phase       = single(niiPhase.img); clear niiPhase;
max_phase   = max(phase(:));
min_phase   = min(phase(:));
phase       = (phase - min_phase) / (max_phase - min_phase) *2*pi - pi;

end