%% output = get_variable_from_headerAndExtraData(headerAndExtraData, variableName)
%
% Input
% --------------
% headerAndExtraData : structure contains any formats of the target variable
% variableName       : name of the variable, string
% matrixSize         : validating matrix size
%
% Output
% --------------
% output             : target variable
%
% Description: get image either directly from structure or specified filename
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 13 August 2021 (v1.0)
% Date modified:
%
%
function output = get_variable_from_headerAndExtraData(headerAndExtraData, variableName, matrixSize)


if isempty(headerAndExtraData.(variableName)) && ~isempty(headerAndExtraData.availableFileList.(variableName))
    % if variable empty but its filename is available
    output = double(load_nii_img_only(headerAndExtraData.availableFileList.(variableName)));
    
    if nargin == 3
        % in case of odd number matrix input
        if ~isequal(size(output), matrixSize(:).')
            output = double(zeropad_odd_dimension(output,'pre')); 
        end
    end
        

else
    % if no fieldmapSD variable provided in any formats
    output = headerAndExtraData.(variableName);
    
end

end