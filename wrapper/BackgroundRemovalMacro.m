%% RDF = BackgroundRemovalMacro(totalField,mask,matrixSize,voxelSize,algorParam,headerAndExtraData)
%
% Input
% --------------
% totalField    : total field map (background + tissue fields), in Hz
% mask          : signal mask
% matrixSize    : size of the input image
% voxelSize     : spatial resolution of each dimension of the data, in mm
% algorParam    : structure contains fields with algorithm-specific parameter(s)
% headerAndExtraData : structure contains extra header info/data for the algorithm
%
% Output
% --------------
% RDF           : local field map
%
% Description: This is a wrapper function to access individual background field removal
%              algorithms for SEPIA (default: 'VSHARP')
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 28 June 2017
% Date modified: 29 September 2017
% Date modified: 1 April 2019
% Date modified: 24 May 2019
% Date modified: 8 March 2020 (v0.8.0)
%
function RDF = BackgroundRemovalMacro(totalField,mask,matrixSize,voxelSize,algorParam,headerAndExtraData)

sepia_universal_variables;
methodBFRName = lower(methodBFRName);

voxelSize       = double(voxelSize(:).');
matrixSize      = double(matrixSize(:).');

algorParam  	= check_and_set_SEPIA_algorithm_default(algorParam);
method          = algorParam.bfr.method;
refine          = algorParam.bfr.refine;
erode_radius	= algorParam.bfr.erode_radius;

headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);

disp('-----------------------------');
disp('Background field removal step');
disp('-----------------------------');

%% zero padding for odd number dimension
fprintf('Zero-padding data if the input images have odd number matrix size...');
totalField  = double(zeropad_odd_dimension(totalField,'pre'));
mask        = double(zeropad_odd_dimension(mask,'pre'));
% additional input
if ~isempty(headerAndExtraData.N_std)
    headerAndExtraData.N_std = double(zeropad_odd_dimension(headerAndExtraData.N_std,'pre'));
end
matrixSize_new = size(totalField);

fprintf('Done!\n');

%% core of background field removal
disp('Removing background field...');
disp(['The following method is being used: ' method]);

for k = 1:length(wrapper_BFR_function)
    if strcmpi(method,methodBFRName{k})
        RDF = feval(wrapper_BFR_function{k},totalField,mask,matrixSize_new,voxelSize,algorParam, headerAndExtraData);
    end
end
disp('Done!');

%% If refine is needed, do it now
if refine
    fprintf('Performing polynomial fitting...');
    % PolyFit required data to be double type
    [~,RDF,~]   = PolyFit(double(RDF),RDF~=0,4);
    fprintf('Done!\n')
end

% get non-zero mask
if erode_radius > 0
    fprintf(['Eroding ' num2str(erode_radius) ' voxel(s) from edges...']);
    maskFinal = RDF ~=0;
    maskFinal = imfill(maskFinal,'holes');
    maskFinal = imerode(maskFinal,strel('sphere',erode_radius));
    % also remove the mask on the edges
    maskFinal(:,:,end-erode_radius:end) = 0;
    maskFinal(:,:,1:erode_radius)       = 0;
    maskFinal(:,end-erode_radius:end,:) = 0;
    maskFinal(:,1:erode_radius,:)       = 0;
    maskFinal(end-erode_radius:end,:,:) = 0;
    maskFinal(1:erode_radius,:,:)       = 0;
    RDF = RDF .* double(maskFinal);
    fprintf('Done!\n')
end

% remove zero padding 
RDF = double(zeropad_odd_dimension(RDF,'post',matrixSize));


end
