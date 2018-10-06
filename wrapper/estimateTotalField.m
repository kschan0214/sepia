%% function [totalField,N_std,fieldmapUnwrapAllEchoes] = estimateTotalField(fieldMap,magn,matrixSize,voxelSize,varargin)
%
% Example:
%   [totalField,N_std] = estimateTotalField(fieldMap,magn,matrixSize,voxelSize,...
%                           'Method','Optimum weights','Unwrap','bestpath3d',...
%                           'unit','ppm','TE',TE,'B0',3,...
%                           'mask',mask);
%   [totalField,N_std] = estimateTotalField(fieldMap,magn,matrixSize,voxelSize,...
%                           'Method','MEDI nonlinear fit','Unwrap','Laplacian','unit','ppm');
%   [totalField,N_std] = estimateTotalField(fieldMap,magn,matrixSize,voxelSize,...
%                           'Method','MEDI nonlinear fit','Unwrap','rg','unit','ppm');
%   [totalField,N_std] = estimateTotalField(fieldMap,magn,matrixSize,voxelSize,...
%                           'Method','MEDI nonlinear fit','Unwrap','gc','unit','ppm','Subsampling',2);
%
% Input
% ----------------
%   fieldMap            : wrapped field map across echoes
%   magn                : magnitude of multi-echo images
%   matrixSize          : images matrix size
%   voxelSize           : images voxel size
%   Name/Value pairs:
%       'method'        -   Method of echo phase combination ('Optimum weights','MEDI nonlinear fit')
%       'Unwarp'        -   phase unwrapping method (Laplacian,rg,gc,bestpath3d)
%       'TE'            -   Echo times
%       'B0'            -   Magnetic field strength
%       'unit'          -   unit of the total field
%       'mask'          -   mask image
%
% Output
% ----------------
%   total field         : unwrapped total field
%   N_std               : noise standard deviation in field map
%
% Description: compute the unwrapped total field from complex-valued
%              multi-echo data
%
% This code is modified from the T2starAndFieldCalc.m from Jose P. Marques
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 31 May, 2017
% Date modified: 27 February 2018
% Date modified: 16 September 2018
%
function [totalField,N_std,fieldmapUnwrapAllEchoes] = estimateTotalField(fieldMap,magn,matrixSize,voxelSize,varargin)
% Larmor frequency of 1H
gyro = 42.57747892;     % (MHz/T)

matrixSize = matrixSize(:).';
voxelSize = voxelSize(:).';
fieldmapUnwrapAllEchoes = [];

%% Parsing argument input flags
if ~isempty(varargin)
    [echoCombine, TE, fieldStrength, unit, unwrapMethod, subsampling, mask] = parse_varargin_Main(varargin);
    if isempty(mask)
        mask = magn(:,:,:,1)>0;
    end
else
    % predefine paramater: if no varargin, use Laplacian
    disp('No method selected. Using the default setting...');
    echoCombine = 'Optimum weights';
    unwrapMethod = 'Laplacian';
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
switch echoCombine
    case 'Optimum weights'
        if size(fieldMap,4) > 1
        %%%%%%%%%%%%%%%%%%%%%%%% Multi-echo %%%%%%%%%%%%%%%%%%%%%%%%
            % compute wrapped phase shift between successive echoes
            fieldMapEchoTemp=angle(exp(1i*fieldMap(:,:,:,2:end))./exp(1i*fieldMap(:,:,:,1:end-1)));
            % unwrap each echo phase shift
            for k=1:size(fieldMapEchoTemp,4)
                tmp = UnwrapPhaseMacro(fieldMapEchoTemp(:,:,:,k),matrixSize,voxelSize,...
                        'method',unwrapMethod,'Magn',magn,'subsampling',subsampling,'mask',mask);
                tmp2(:,:,:,k)=tmp-round(tmp(pos(1),pos(2),pos(3))/(2*pi))*2*pi;
            end
            % get phase accumulation over all echoes
            phaseShiftUnwrapAllEchoes=cumsum(tmp2,4);    

            % compute unwrapped phase shift between successive echoes
            fieldMapUnwrap = zeros(size(phaseShiftUnwrapAllEchoes));
            for k=1:size(phaseShiftUnwrapAllEchoes,4)
                fieldMapUnwrap(:,:,:,k)=phaseShiftUnwrapAllEchoes(:,:,:,k)/(TE(k+1)-TE(1));
            end

            % Robinson et al. 2017 NMR Biomed Appendix A2
            fieldMapSD = zeros(size(phaseShiftUnwrapAllEchoes));
            for k=1:size(phaseShiftUnwrapAllEchoes,4)
                fieldMapSD(:,:,:,k) = 1./(TE(k+1)-TE(1)) ...
                    * sqrt((magn(:,:,:,1).^2+magn(:,:,:,k+1).^2)./((magn(:,:,:,1).*magn(:,:,:,k+1)).^2));
            end
            % weights are inverse of the field map variance
            weight = bsxfun(@rdivide,1/(fieldMapSD.^2),sum(1/(fieldMapSD.^2),4));
            % Weighted average of unwrapped phase shift
            totalField = sum(fieldMapUnwrap .* weight,4);
            totalField(isnan(totalField))=0;
            totalField(isinf(totalField))=0;
            % standard deviation of field map from weighted avearging
            totalFieldVariance = sum(weight.^2 .* fieldMapSD.^2,4);
            totalFieldSD = sqrt(totalFieldVariance);
            totalFieldSD(isnan(totalFieldSD)) = 0;
            totalFieldSD(isinf(totalFieldSD)) = 0;
            totalFieldSD = totalFieldSD./norm(totalFieldSD(totalFieldSD~=0));
            
            % get the unwrapped phase accumulation across echoes
            % unwrap first echo
            tmp = UnwrapPhaseMacro(fieldMap(:,:,:,1),matrixSize,voxelSize,...
                'method',unwrapMethod,'Magn',magn,'subsampling',subsampling,'mask',mask);
            tmp = tmp-round(tmp(pos(1),pos(2),pos(3))/(2*pi))*2*pi;
            fieldmapUnwrapAllEchoes = cat(4,tmp,tmp2);
            fieldmapUnwrapAllEchoes = cumsum(fieldmapUnwrapAllEchoes,4);
        else
        %%%%%%%%%%%%%%%%%%%%%%%% Single echo %%%%%%%%%%%%%%%%%%%%%%%%
            tmp = UnwrapPhaseMacro(fieldMap,matrixSize,voxelSize,...
                        'method',unwrapMethod,'Magn',magn,'subsampling',subsampling,'mask',mask);
            tmp = tmp-round(tmp(pos(1),pos(2),pos(3))/(2*pi))*2*pi;
            totalField = tmp/(TE(1));
            totalFieldSD = 1./magn;
            totalFieldSD(isnan(totalFieldSD)) = 0;
            totalFieldSD(isinf(totalFieldSD)) = 0;
            totalFieldSD = totalFieldSD./norm(totalFieldSD(totalFieldSD~=0));
        end
        
        N_std = totalFieldSD;
        
    case 'MEDI nonlinear fit'
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
                        'method',unwrapMethod,'Magn',rss_magn,'subsampling',subsampling,'mask',mask) / dt;
end

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
function [echoCombine,TE,fieldStrength,unit,unwrapMethod,subsampling,mask] = parse_varargin_Main(arg)
echoCombine = 'Optimum weights';
TE = 1;
fieldStrength = 3;
subsampling = 1;
unit = 'ppm';
unwrapMethod = 'laplacian';
mask = [];
for kkvar = 1:length(arg)
    if strcmpi(arg{kkvar},'method')
        echoCombine = arg{kkvar+1};
        continue
    end
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
        unwrapMethod = lower(arg{kkvar+1});
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
