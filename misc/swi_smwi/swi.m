%% [pSWI,dSWI,swi_phase] = swi(magn,phase,filterSize,thres,m)
%
% Input
% --------------
% magn          : image that combines with phase contrast, can be 3D/4D
% phase         : phase map (in rad)
% filterSize    : window length of Hamming filter
% thres         : threshold used to create positive/negative phase contrast
% m             : power to enhance phase contrast
% method        : method of process phase map ('default' or 'multiecho')
%
% Output
% --------------
% pSWI          : positive phase enahcned SWI
% nSWI          : negative phase enahcned SWI
% swi_phase     : high-pass filtered phase
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 14 April 2019
% Date last modified:
%
%
function [pSWI,nSWI,swi_phase] = swi(magn,phase,filterSize,thres,m,method)
% check and set default
if nargin < 6 || size(magn,4) == 1
    method = 'default';
end
if nargin < 5
    m = 4;
end
if nargin < 4 || thres < 0
    warning('Threshold has to be greater than zero. Using pi instead...');
    thres = pi;
end
% convert input to double
magn    = double(magn);
phase   = double(phase);
m       = round(m);

switch method
    case 'default'
        
    case 'multiecho'
        phaseDiff = angle(exp(1i*phase(:,:,:,2:end))./exp(1i*phase(:,:,:,1:end-1)));
        phaseDiff = cat(4,phase(:,:,:,1),phaseDiff);
        phase = phaseDiff;
end

% display algorithm parameters
disp(['Filter size: ' num2str(filterSize)]);
disp(['Threshold (rad): ' num2str(thres)]);
disp(['Contrast power: ' num2str(m)]);
disp(['Method to process phase map: ' method]);

% convert phase to complex value
phase_cplx = exp(1i*phase);

dim = size(phase_cplx);

% for single echo
if length(dim) < 4
    dim(4) = 1;
end

% create a 2D low-pass filtre
lowpass_filter = hamming(filterSize)*hamming(filterSize)'; 

swi_phase = zeros(dim);
% loops all echoes
for kt = 1:dim(4)
    disp(['Processing echo ' num2str(kt) ' ...']);
    % loops all slices, the thresholding is actually performed on
    % slice-by-slice basis
    for kz = 1:dim(3)
        % get slice with complex-valued data
        c = phase_cplx(:,:,kz, kt);
        
        % filtre the complex-valued slice with the above low-pass filtre
        c_lowpass = filter2(lowpass_filter, c);   
        
        % complex-valued division between the original data and low-passed
        % data is equivalent to high-pass filtred the original data
        swi_highpass = c ./ c_lowpass;
        % KC: compute the high-pass filtred phase from the complex-valued slice
        swi_phase(:, :, kz, kt) = angle(swi_highpass);
    end 
end

switch method
    case 'default'
        
    case 'multiecho'
        swi_phase = cumsum(swi_phase,4);
end

% positive phase mask
pSmask                  = (thres-swi_phase)/thres;
pSmask(swi_phase>thres)	= 0; 
pSmask(swi_phase<0)   	= 1;
% negative phase mask
dSmask                 	= (thres+swi_phase)/thres;
dSmask(swi_phase>0)   	= 1; 
dSmask(swi_phase<-thres)= 0;

pSWI = bsxfun(@times,magn,pSmask.^m);
nSWI = bsxfun(@times,magn,dSmask.^m);

end