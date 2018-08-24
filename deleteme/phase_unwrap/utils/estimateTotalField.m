%% function [totalField,N_std] = estimateTotalField(fieldMap,magn,matrixSize,voxelSize,varargin)
%
% Usage:
%   [totalField,N_std] = estimateTotalField(fieldMap,magn,matrixSize,voxelSize,...
%                           'Unwrap','Laplacian','unit','ppm');
%   [totalField,N_std] = estimateTotalField(fieldMap,magn,matrixSize,voxelSize,...
%                           'Unwrap','rg','unit','ppm');
%   [totalField,N_std] = estimateTotalField(fieldMap,magn,matrixSize,voxelSize,...
%                           'Unwrap','gc','unit','ppm','Subsampling',2);
%   [totalField,N_std] = estimateTotalField(fieldMap,magn,matrixSize,voxelSize,...
%                           'Unwrap','bestpath3d','unit','ppm','TE',TE,'B0',3,...
%                           'mask',mask);
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
% This code is modified from the T2starAndFieldCalc.m from Jose P. Marques
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 31 May, 2017
% Date modified: 27 February 2018
% Date modified: 27 May 2018
%
function [totalField,N_std] = estimateTotalField(fieldMap,magn,matrixSize,voxelSize,varargin)
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
end

% find the centre of mass
pos=round(centerofmass(magn));

%% Core
if size(fieldMap,4) > 1
    % Multi-echo
    fieldMapEchoTemp=angle(exp(1i*fieldMap(:,:,:,2:end))./exp(1i*fieldMap(:,:,:,1:end-1)));
    for k=1:size(fieldMapEchoTemp,4)
        tmp = UnwrapPhaseMacro(fieldMapEchoTemp(:,:,:,k),matrixSize,voxelSize,...
                'method',method,'Magn',magn,'subsampling',subsampling,'mask',mask);
        tmp2(:,:,:,k)=tmp-round(tmp(pos(1),pos(2),pos(3))/(2*pi))*2*pi;
    end
    fieldMapTemp=cumsum(tmp2,4);
    % fieldMapTemp=gather(fieldMapTemp);
    for k=1:size(fieldMapTemp,4)
        fieldMapTemp(:,:,:,k)=fieldMapTemp(:,:,:,k)/(TE(k+1)-TE(1));
    end
    % denominator=zeros(matrixSize(1:3));
    % numerator=zeros(matrixSize(1:3));
    fieldMapSD = zeros(size(fieldMapTemp));

    for k=1:size(fieldMapTemp,4)
    %     weight_k(:,:,:,k) = 1./((sqrt((magn(:,:,:,1).^2+magn(:,:,:,k+1).^2)./((magn(:,:,:,1).*magn(:,:,:,k+1)).^2))/(abs(TE(k+1)-TE(1)))).^2);
        fieldMapSD(:,:,:,k) = 1./(TE(k+1)-TE(1)) ...
            * sqrt((magn(:,:,:,1).^2+magn(:,:,:,k+1).^2)./((magn(:,:,:,1).*magn(:,:,:,k+1)).^2));
    %     numerator = numerator + fieldMapTemp(:,:,:,k).*weight_k;
    %     denominator = denominator + weight_k;
    %     numerator= numerator + (fieldMapTemp(:,:,:,k)) ...
    %         .*abs((TE(k+1)-TE(1))*magn(:,:,:,k+1)).^2 ...
    %         ./(abs(magn(:,:,:,k+1)).^2+abs(magn(:,:,:,1)).^2);
    %     denominator= denominator + ...
    %         abs((TE(k+1)-TE(1))*magn(:,:,:,k+1)).^2 ...
    %         ./(abs(magn(:,:,:,k+1)).^2+abs(magn(:,:,:,1)).^2);
    end
    weight = bsxfun(@rdivide,1/(fieldMapSD.^2),sum(1/(fieldMapSD.^2),4));
    totalField = sum(fieldMapTemp .* weight,4);
    % totalField = sum(fieldMapTemp .* (1./fieldMapSD.^2),4)./sum(1./fieldMapSD.^2,4);
    % totalField = numerator./denominator;
    totalField(isnan(totalField))=0;
    totalField(isinf(totalField))=0;
    % standard deviation of weighted mean
    % totalFieldSD = sqrt(1./sum(1./fieldMapSD.^2,4));
    totalFieldVariance = sum(weight.^2 .* fieldMapSD.^2,4);
    totalFieldSD = sqrt(totalFieldVariance);
    totalFieldSD(isnan(totalFieldSD)) = 0;
    totalFieldSD(isinf(totalFieldSD)) = 0;
    totalFieldSD = totalFieldSD./norm(totalFieldSD(totalFieldSD~=0));
else
    % Single echo
    tmp = UnwrapPhaseMacro(fieldMap,matrixSize,voxelSize,...
                'method',method,'Magn',magn,'subsampling',subsampling,'mask',mask);
    tmp = tmp-round(tmp(pos(1),pos(2),pos(3))/(2*pi))*2*pi;
    totalField = tmp/(TE(1));
    totalFieldSD = 1./magn;
    totalFieldSD(isnan(totalFieldSD)) = 0;
    totalFieldSD(isinf(totalFieldSD)) = 0;
    totalFieldSD = totalFieldSD./norm(totalFieldSD(totalFieldSD~=0));
end

% totalField now in rads^-1, matching tmp2 to the smae unit
% tmp2 = tmp2/dTE;

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

N_std = totalFieldSD;
% N_std = sqrt(mean(bsxfun(@minus,tmp2,totalField).^2,4));

end

%% find the centre of mass
function coord=centerofmass(data)
data=abs(data);
dims=size(data);
    for k=1:length(dims)
    %     datatemp=permute(data,[k ]);
    dimsvect=ones([1, length(dims)]);
    dimsvect(k)=dims(k);
    temp=bsxfun(@times,(data),reshape(1:dims(k),dimsvect));
    coord(k)=sum(temp(:))./sum(data(:));
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
