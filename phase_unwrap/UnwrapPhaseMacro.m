%% function unwrappedField = UnwrapPhaseMacro(wrappedField,matrixSize,voxelSize,varargin)
%
% Usage: 
%   unwrappedField = UnwrapPhaseMacro(wrappedField,matrixSize,voxelSize,...
%                       'method','laplacian');
%   unwrappedField = UnwrapPhaseMacro(wrappedField,matrixSize,voxelSize,...
%                       'method','rg','Magn',magn);
%   unwrappedField = UnwrapPhaseMacro(wrappedField,matrixSize,voxelSize,...
%                       'method','gc','Magn',magn,'subsampling',2);
%   unwrappedField = UnwrapPhaseMacro(wrappedField,matrixSize,voxelSize,...
%                       'method','bestpath3d','mask',mask);
%
% Description: Wrapper for phase unwrapping (default using Laplacian)
%   Flags:
%       'method'        : phase unwrapping method, 
%                          'laplacian', 'laplacian_stisuite','rg', 'gc' and
%                          'bestpath3d' 
%
%       Laplacian
%       ----------------
%
%       rg (RegionGrowing)
%       ----------------           
%       'Magn'          : magnitude data
%
%       gc (Graphcut)
%       ----------------
%       'Magn'          : magnitude data
%       'subsampling'	: Downsampling factor for speed
%
%       bestpath3d (3D best path)
%       ----------------
%       'mask'          : brain mask
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 29 June 2017
% Date last modified: 24 May 2018
%
function unwrappedField = UnwrapPhaseMacro(wrappedField,matrixSize,voxelSize,varargin)

matrixSize = matrixSize(:).';
voxelSize = voxelSize(:).';

%% Parsing argument input flags
if ~isempty(varargin)
    for kvar = 1:length(varargin)
        if strcmpi(varargin{kvar},'method')
            switch lower(varargin{kvar+1})
                case 'laplacian'
                    method = 'Laplacian';
                    break
                case 'laplacian_stisuite'
                    method = 'Laplacian_stisuite';
                    break
                case 'rg'
                    method = 'RegionGrowing';
                    [magn] = parse_varargin_RegionGrowing(varargin);
                    if isempty(magn)
                        disp('Running algorithm without magnitude image could be problematic');
                        magn = ones(matrixSize);
                    end
                    break
                case 'gc'
                    method = 'Graphcut';
                    [magn, subsampling] = parse_varargin_Graphcut(varargin);
                    if isempty(magn)
                        disp('Running algorithm without magnitude image could be problematic');
                        magn = ones(matrixSize);
                    end
                    break
                case 'bestpath3d'
                    method = 'BestPath3D';
                    [mask] = parse_varargin_UnwrapPhase_3DBestPath(varargin);
                    if isempty(mask)
                        disp('Running algorithm without brain mask could be problematic');
                        mask = ones(matrixSize);
                    end
                    break
            end
        end
    end
else
    % predefine paramater: if no varargin, use Laplacian
    disp('No method selected. Using default setting.');
    method = 'Laplacian';
end

disp(['The following unwrapping method is being used: ' method]);
%% phase unwrapping
switch method
    case 'Laplacian'
        unwrappedField = unwrapLaplacian(wrappedField,matrixSize,voxelSize);
    case 'Laplacian_stisuite'
        unwrappedField = MRPhaseUnwrap(wrappedField,'voxelsize',voxelSize,'padsize',[12,12,12]);
    case 'RegionGrowing'
        unwrappedField = unwrapPhase(magn,wrappedField,matrixSize);
    case 'Graphcut'
        disp(['Graphcut subsampling factor: ' num2str(subsampling)]);
        unwrappedField = unwrapping_gc(wrappedField,magn,voxelSize,subsampling);
    case 'BestPath3D'
        try
            unwrappedField = UnwrapPhase_3DBestPath(wrappedField,mask,matrixSize);
        catch
            disp('The library cannot be run in this platform, running Laplacian unwrapping instead...');
            unwrappedField = unwrapLaplacian(wrappedField,matrixSize,voxelSize);
        end
end

end
