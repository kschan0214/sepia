%% function unwrappedField = UnwrapPhaseMacro(wrappedField,matrixSize,voxelSize,varargin)
%
% Usage: 
%   unwrappedField = UnwrapPhaseMacro(wrappedField,matrixSize,voxelSize,...
%                       'method','laplacian');
%   unwrappedField = UnwrapPhaseMacro(wrappedField,matrixSize,voxelSize,...
%                       'method','regiongrowing','Magn',magn);
%   unwrappedField = UnwrapPhaseMacro(wrappedField,matrixSize,voxelSize,...
%                       'method','graphcut','Magn',magn,'subsampling',2);
%   unwrappedField = UnwrapPhaseMacro(wrappedField,matrixSize,voxelSize,...
%                       'method','graphcut','mask',mask);
%
% Description: Wrapper for phase unwrapping (default using Laplacian)
%   Flags:
%       'method'        : phase unwrapping method, 
%                          'Laplacian', 'RegionGrowing' and 'Graphcut' 
%
%       Laplacian
%       ----------------
%
%       RegionGrowing
%       ----------------           
%       'Magn'          : magnitude data
%
%       Graphcut
%       ----------------
%       'Magn'          : magnitude data
%       'subsampling'	: Downsampling factor for speed
%
%       Jena
%       ----------------
%       'mask'          : brain mask
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 29 June 2017
% Date last modified:
%
function unwrappedField = UnwrapPhaseMacro(wrappedField,matrixSize,voxelSize,varargin)
%% Parsing argument input flags
if ~isempty(varargin)
    for kvar = 1:length(varargin)
        if strcmpi(varargin{kvar},'method')
            switch lower(varargin{kvar+1})
                case 'laplacian'
                    method = 'Laplacian';
                    break
                case 'regiongrowing'
                    [method, magn] = parse_vararginRegionGrowing(varargin);
                    if isempty(magn)
                        magn = ones(matrixSize);
                    end
                    break
                case 'graphcut'
                    [method, magn, subsampling] = parse_vararginGraphcut(varargin);
                    if isempty(magn)
                        magn = ones(matrixSize);
                    end
                    break
                case 'jena'
                    [method, mask] = parse_vararginJena(varargin);
                    if isempty(mask)
                        mask = ones(matrixSize);
                    end
                    break
            end
        end
    end
else
    % predefine paramater: if no varargin, use Laplacian
    disp('No method selected. Using the default setting:');
    method = 'Laplacian'
end

%% phase unwrapping
switch method
    case 'Laplacian'
        unwrappedField = unwrapLaplacian(wrappedField,matrixSize,voxelSize);
    case 'RegionGrowing'
        unwrappedField = unwrapPhase(magn,wrappedField,matrixSize);
    case 'Graphcut'
        unwrappedField = unwrapping_gc(wrappedField,magn,voxelSize,subsampling);
    case 'Jena'
        unwrappedField = unwrapJena(wrappedField,mask,matrixSize);
end

end

%% Parsing varargin
% Region growing
function [method, magn] = parse_vararginRegionGrowing(arg)
method = 'RegionGrowing';
magn = [];
for kkvar = 1:length(arg)
    if strcmpi(arg{kkvar},'Magn')
        magn = arg{kkvar+1};
        continue
    end
end
end

% Graphcut
function [method, magn, subsampling] = parse_vararginGraphcut(arg)
method = 'Graphcut';
magn = [];
subsampling = 1;
for kkvar = 1:length(arg)
    if strcmpi(arg{kkvar},'Magn')
        magn = arg{kkvar+1};
        continue
    end
    if  strcmpi(arg{kkvar},'Subsampling')
        subsampling = arg{kkvar+1};
        continue
    end
end
end

% Jena
function [method, mask] = parse_vararginJena(arg)
method = 'Jena';
mask = [];
for kkvar = 1:length(arg)
    if strcmpi(arg{kkvar},'mask')
        mask = arg{kkvar+1};
        continue
    end
end
end
