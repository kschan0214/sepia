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
% Date modified: 27 May 2018
% Date modified: 24 May 2019
% Date modified: 6 March 2020 (v0.8.0)
%
function unwrappedField = UnwrapPhaseMacro(wrappedField,matrixSize,voxelSize,varargin)

sepia_universal_variables;
methodUnwrapName = lower(methodUnwrapName);

matrixSize  = matrixSize(:).';
voxelSize   = voxelSize(:).';
mask        = [];

%% Parsing argument input flags
if ~isempty(varargin)
    for kvar = 1:length(varargin)
        if strcmpi(varargin{kvar},'method')
            method = varargin{kvar+1};
            switch lower(method)
                case methodUnwrapName{1}    % laplacian medi

                case methodUnwrapName{2}    % laplacian sti suite
                    
                case methodUnwrapName{3}    % 3d best path

                case methodUnwrapName{4}    % region growing medi
                    
                    [magn] = parse_varargin_RegionGrowing(varargin);
                    if isempty(magn)
                        disp('Running algorithm without magnitude image could be problematic');
                        magn = ones(matrixSize);
                    end

                case methodUnwrapName{5}    % graphcut
                    
                    [magn, subsampling] = parse_varargin_Graphcut(varargin);
                    if isempty(magn)
                        disp('Running algorithm without magnitude image could be problematic');
                        magn = ones(matrixSize);
                    end
                    
                case methodUnwrapName{6}    % segue

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

% give a warning if no mask is provided
if isempty(mask)
    mask = ones(matrixSize);
    warning('Running algorithm without brain mask could be problematic');
end

% use same data type (20190529: single-precision can affect the result of
% some methods)
wrappedField = double(wrappedField);
mask         = double(mask);
matrixSize   = double(matrixSize);
voxelSize    = double(voxelSize);
if exist('magn','var')
    magn = double(magn);
end

disp(['The following unwrapping method is being used: ' method]);
%% phase unwrapping
% add method path
sepia_addpath(method);

% Laplacian based method required prior zero padding for odd number dimension
if strcmpi(method,methodUnwrapName{1}) || strcmpi(method,methodUnwrapName{2})
    wrappedField = zeropad_odd_dimension(wrappedField,'pre');
    matrixSize_new = size(wrappedField);
end

switch method
    case methodUnwrapName{1}    % laplacian medi
        
        % Laplacian unwrapping
        unwrappedField = unwrapLaplacian(wrappedField,matrixSize_new,voxelSize);
        
    case methodUnwrapName{2}    % laplacian sti suite

        unwrappedField = MRPhaseUnwrap(wrappedField,'voxelsize',voxelSize,'padsize',[12,12,12]);
        
    case methodUnwrapName{4}    % region growing
        if size(magn,4) > 1
            magn = sqrt(sum(abs(magn).^2,4));
        end
        magn = magn .* mask;
        unwrappedField = unwrapPhase(magn,wrappedField,matrixSize);
    case methodUnwrapName{5}    % graphcut
        disp(['Graphcut subsampling factor: ' num2str(subsampling)]);
        if size(magn,4) > 1
            magn = sqrt(sum(abs(magn).^2,4));
        end
        magn = magn .* mask;
        unwrappedField = unwrapping_gc(wrappedField,magn,voxelSize,subsampling);
    case methodUnwrapName{3}    % 3d best path
        try
            unwrappedField = UnwrapPhase_3DBestPath(wrappedField,mask,matrixSize);
            
        catch ME
            sepia_addpath(methodUnwrapName{4});
            
            warning('The library cannot be run in this platform, running region growing unwrapping instead...');
            [magn] = parse_varargin_RegionGrowing(varargin);
            if isempty(magn)
                disp('Running algorithm without magnitude image could be problematic');
                magn = ones(matrixSize,'like',matrixSize);
            end
            unwrappedField = unwrapPhase(magn,wrappedField,matrixSize);
        end
    case methodUnwrapName{6}    % segue
        try 
            Inputs.Mask     = mask;
            Inputs.Phase    = wrappedField;
            unwrappedField  = SEGUE(Inputs);
        catch ME
            sepia_addpath(methodUnwrapName{4});
            
            warning('Problem suing function. Running MEDI region growing unwrapping instead...');
            [magn] = parse_varargin_RegionGrowing(varargin);
            if isempty(magn)
                disp('Running algorithm without magnitude image could be problematic');
                magn = ones(matrixSize,'like',matrixSize);
            end
            unwrappedField = unwrapPhase(magn,wrappedField,matrixSize);
        end
end

% remove zero padding with Laplacian based method result
if strcmpi(method,methodUnwrapName{1}) || strcmpi(method,methodUnwrapName{2})
    unwrappedField = zeropad_odd_dimension(unwrappedField,'post',matrixSize);
end

unwrappedField = double(unwrappedField);

end
