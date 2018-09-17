%% function [totalField,N_std] = estimateTotalField_MEDI(fieldMap,magn,matrixSize,voxelSize,varargin)
%
% Description: compute the unwrapped total field from complex-valued
%              multi-echo data
%
% Input
% ----------------
%   fieldMap            : wrapped field map across echoes
%   magn                : magnitude of multi-echo images
%   matrixSize          : images matrix size
%   voxelSize           : images voxel size
%   Flags:
%       'TE'            -   Echo times
%       'B0'            -   Magnetic field strength
%       'unit'          -   unit of the total field
%       'Unwarp'        -   phase unwrapping method (Laplacian,rg,gc,bestpath3d)
%
% Output
% ----------------
%   total field         : unwrapped total field
%   N_std               : noise standard deviation in field map
%
% This code is modified from MEDI toolbox README.m
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 27 May 2018
% Date modified:  
%
function [totalField,N_std] = estimateTotalField_MEDI(fieldMap,magn,matrixSize,voxelSize,varargin)
% Larmor frequency of 1H
gyro = 42.57747892;     % (MHz/T)

matrixSize = matrixSize(:).';
voxelSize = voxelSize(:).';

%% Parsing argument input flags
if ~isempty(varargin)
    [TE, fieldStrength, unit, method, subsampling, mask] = parse_varargin_Main(varargin);
    if isempty(mask)
        mask = magn(:,:,:,1)>0;
    end
else
    % predefine paramater: if no varargin, use Laplacian
    disp('No method selected. Using the default setting...');
    method = 'Laplacian';
    TE = 1;
    fieldStrength = 3;
    unit = 'ppm';
    subsampling = 1;
    mask = magn(:,:,:,1)>0;
end

if length(TE)>1
    dt = TE(2)-TE(1);
else
    dt = TE;
end

%% Core
sepia_addpath('nonlinearfit');
if numel(TE)>1 && ((TE(2)-TE(1))-(TE(3)-TE(2))>1e-5)
    % Estimate the frequency offset in each of the voxel using a complex
    % fitting (uneven echo spacing)
    [iFreq_raw, N_std] = Fit_ppm_complex_TE(magn.*exp(-1i*fieldMap),TE);
else
    % Estimate the frequency offset in each of the voxel using a complex
    % fitting (even echo spacing)
    [iFreq_raw, N_std] = Fit_ppm_complex(magn.*exp(-1i*fieldMap));
end

% Compute magnitude image
rss_magn = sqrt(sum(abs(magn).^2,4));

% Spatial phase unwrapping
totalField = UnwrapPhaseMacro(iFreq_raw,matrixSize,voxelSize,...
                'method',method,'Magn',rss_magn,'subsampling',subsampling,'mask',mask) / dt;

disp(['The resulting field map with the following unit: ' unit]);
switch unit
    case 'ppm'
        totalField = (totalField/(2*pi))/(fieldStrength*gyro);
%         tmp2 = (tmp2/(2*pi))/(fieldStrength*gamma);
    case 'rad'
        totalField = totalField*dt;
%         tmp2 = tmp2*dt;
    case 'hz'
        totalField = totalField/(2*pi);
%         tmp2 = tmp2/(2*pi);
    case 'radhz'
    otherwise
        disp(['Input unit is invalid. radHz is used instead.']);
end

end

%% Parsing varargin
function [TE,fieldStrength,unit,method,subsampling,mask] = parse_varargin_Main(arg)
TE = 1;
fieldStrength = 3;
subsampling = 1;
unit = 'ppm';
method = 'laplacian';
mask = [];
for kkvar = 1:length(arg)
    if strcmpi(arg{kkvar},'TE')
        TE = arg{kkvar+1};
        continue
    end
    if  strcmpi(arg{kkvar},'B0')
        fieldStrength = arg{kkvar+1};
        continue
    end
    if  strcmpi(arg{kkvar},'unit')
        unit = lower(arg{kkvar+1});
        continue
    end
    if  strcmpi(arg{kkvar},'Unwrap')
        method = lower(arg{kkvar+1});
        continue
    end
    if  strcmpi(arg{kkvar},'Subsampling')
        subsampling = arg{kkvar+1};
        continue
    end
    if  strcmpi(arg{kkvar},'mask')
        mask = arg{kkvar+1};
        continue
    end
end
end
