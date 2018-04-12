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
% Date last modified:
%
%
function phase = DICOM2Phase(dicomPhase)
% Siemens DICOM phase range from -4069 to 4094
% Scale the full range to [0,1]
dicomPhase = (dicomPhase +4096)/(4094+4096);
% Scale to [0,2*pi]
dicomPhase = dicomPhase * 2*pi;
% shift to [-pi,pi)
phase = dicomPhase -pi;
end