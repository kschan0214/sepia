%% unwrappedField = UnwrapPhaseMacro(wrappedField,mask,matrixSize,voxelSize,algorParam,headerAndExtraData)
%
% Input
% --------------
% wrappedField  : original wrapped phase image, in rad
% mask          : signal mask
% matrixSize    : size of the input image
% voxelSize     : spatial resolution of each dimension of the data, in mm
% algorParam    : structure contains fields with algorithm-specific parameter(s)
% headerAndExtraData : structure contains extra header info/data for the algorithm
%
% Output
% --------------
% unwrappedField : unwrapped phase image
%
% Description: This is a wrapper function to access individual phase
%              unwrapping algorithms for SEPIA (default: 'Laplacian (MEDI)')
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 29 June 2017
% Date modified: 27 May 2018
% Date modified: 24 May 2019
% Date modified: 9 March 2020 (v0.8.0)
% Date modified: 13 August 2021 (v1.0)
%
function unwrappedField = UnwrapPhaseMacro(wrappedField,mask,matrixSize,voxelSize,algorParam,headerAndExtraData)

sepia_universal_variables;
methodUnwrapName = lower(methodUnwrapName);

matrixSize  = double(matrixSize(:).');
voxelSize   = double(voxelSize(:).');

algorParam          = check_and_set_SEPIA_algorithm_default(algorParam);
method              = algorParam.unwrap.unwrapMethod;

headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);

disp('------------------------');
disp('Spatial phase unwrapping');
disp('------------------------');

% give a warning if no mask is provided
if isempty(mask)
    mask = ones(matrixSize);
    warning('Running algorithm without brain mask could be problematic');
end

if isempty(headerAndExtraData.magnitude) && isempty(headerAndExtraData.availableFileList.magnitude) 
    warning('Running algorithm without magnitude image could be problematic');
%     headerAndExtraData.magnitude = ones(matrixSize,'like',matrixSize);
end

%% Laplacian based method required prior zero padding for odd number dimension
% use same data type (20190529: single-precision can affect the result of
% some methods)
fprintf('Zero-padding data if the input images have odd number matrix size...');
if strcmpi(method,methodUnwrapName{1}) || strcmpi(method,methodUnwrapName{2})
    wrappedField    = double(zeropad_odd_dimension(wrappedField,'pre'));
    mask            = double(zeropad_odd_dimension(mask,'pre'));
    % additional input
    if ~isempty(headerAndExtraData.magnitude)
        headerAndExtraData.magnitude = double(zeropad_odd_dimension(headerAndExtraData.magnitude,'pre'));
    end
end
matrixSize_new = size(wrappedField);

fprintf('Done!\n');

%% phase unwrapping
disp('Unwrapping phase image...');
disp(['The following method is being used: ' method]);

for k = 1:length(wrapper_Unwrap_function)
    if strcmpi(method,methodUnwrapName{k})
        try 
            unwrappedField = feval(wrapper_Unwrap_function{k},wrappedField,mask,matrixSize_new,voxelSize,algorParam, headerAndExtraData);
        catch
            % any problem occur run region growing instead
            warning('Problem using this algorithm. Running MEDI region growing unwrapping instead...');
            unwrappedField = feval(wrapper_Unwrap_function{4},wrappedField,mask,matrixSize_new,voxelSize,algorParam, headerAndExtraData);
        end
    end
end
disp('Done!');

%% remove zero padding with Laplacian based method result
if strcmpi(method,methodUnwrapName{1}) || strcmpi(method,methodUnwrapName{2})
    unwrappedField = double(zeropad_odd_dimension(unwrappedField,'post',matrixSize));
end

end
