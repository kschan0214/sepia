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
% Date last modified: 27 May 2018
%
function unwrappedField = UnwrapPhaseMacro(wrappedField,matrixSize,voxelSize,varargin)

matrixSize = matrixSize(:).';
voxelSize = voxelSize(:).';
mask = [];

%% Parsing argument input flags
if ~isempty(varargin)
    for kvar = 1:length(varargin)
        if strcmpi(varargin{kvar},'method')
            switch lower(varargin{kvar+1})
                case 'laplacian'
                    method = 'Laplacian';
%                     break
                case 'laplacian_stisuite'
                    method = 'Laplacian_stisuite';
%                     break
                case 'rg'
                    method = 'RegionGrowing';
                    [magn] = parse_varargin_RegionGrowing(varargin);
                    if isempty(magn)
                        disp('Running algorithm without magnitude image could be problematic');
                        magn = ones(matrixSize);
                    end
%                     break
                case 'gc'
                    method = 'Graphcut';
                    [magn, subsampling] = parse_varargin_Graphcut(varargin);
                    if isempty(magn)
                        disp('Running algorithm without magnitude image could be problematic');
                        magn = ones(matrixSize);
                    end
%                     break
                case 'bestpath3d'
                    method = 'BestPath3D';
%                     [mask] = parse_varargin_UnwrapPhase_3DBestPath(varargin);
%                     if isempty(mask)
%                         disp('Running algorithm without brain mask could be problematic');
%                         mask = ones(matrixSize);
%                     end
%                     break
            end
        end
        if strcmpi(varargin{kvar},'mask')
            mask = varargin{kvar+1};
        end
    end
else
    % predefine paramater: if no varargin, use Laplacian
    disp('No method selected. Using default setting.');
    method = 'Laplacian';
end

if isempty(mask)
    mask = ones(matrixSize);
    warning('Running algorithm without brain mask could be problematic');
end

% add path
sepia_addpath(method);

disp(['The following unwrapping method is being used: ' method]);
%% phase unwrapping
switch method
    case 'Laplacian'
        % check odd matrix dimension
        wrappedField = DataValidation(wrappedField,'pre');
        matrixSize_new = size(wrappedField);
        
        % Laplacian unwrapping
        unwrappedField = unwrapLaplacian(wrappedField,matrixSize_new,voxelSize);
        
        % remove zero-padding if any 
        unwrappedField = DataValidation(unwrappedField,'post',matrixSize);
        
    case 'Laplacian_stisuite'
        % check odd matrix dimension
        wrappedField = DataValidation(wrappedField,'pre');
        
        unwrappedField = MRPhaseUnwrap(wrappedField,'voxelsize',voxelSize,'padsize',[12,12,12]);
        
        % remove zero-padding if any 
        unwrappedField = DataValidation(unwrappedField,'post',matrixSize);
        
    case 'RegionGrowing'
        if size(magn,4) > 1
            magn = sqrt(sum(abs(magn).^2,4));
        end
        magn = magn .* mask;
        unwrappedField = unwrapPhase(magn,wrappedField,matrixSize);
    case 'Graphcut'
        disp(['Graphcut subsampling factor: ' num2str(subsampling)]);
        if size(magn,4) > 1
            magn = sqrt(sum(abs(magn).^2,4));
        end
        magn = magn .* mask;
        unwrappedField = unwrapping_gc(wrappedField,magn,voxelSize,subsampling);
    case 'BestPath3D'
        try
            unwrappedField = UnwrapPhase_3DBestPath(wrappedField,mask,matrixSize);
        catch ME
            warning('The library cannot be run in this platform, running region growing unwrapping instead...');
            [magn] = parse_varargin_RegionGrowing(varargin);
            if isempty(magn)
                disp('Running algorithm without magnitude image could be problematic');
                magn = ones(matrixSize);
            end
            unwrappedField = unwrapPhase(magn,wrappedField,matrixSize);
        end
end

end

%% make sure the size of the input matrix is an even number
function output = DataValidation(input,mode,matrixSize_o)
matrixSize = size(input);

% determine if a dimension needs to be zeropadded
padsize     = zeros(size(matrixSize));
for kd = 1:length(matrixSize)
    if mod(matrixSize(kd),2) == 1
        padsize(kd) = 1;
    end
end

switch mode
    case 'pre'
        % zero padding if the dimension of the matrix is an odd number
        output = padarray(input, padsize, 0,'post');
        
    case 'post'
        % remove zero padding 
        output = input(1:matrixSize_o(1),1:matrixSize_o(2),1:matrixSize_o(3));
        
end
end
