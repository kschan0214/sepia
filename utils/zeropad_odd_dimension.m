%% output = zeropad_odd_dimension(input,mode,matrixSize_o)
%
% Usage:
%
% Input
% --------------
% input         : input image
% mode          : 'pre' - preprocessing (zero padding); 'post' -
%                 postprocessing (remove zero padding)
% matrixSize_o  : original matrix size beofre zero padding
%
% Output
% --------------
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 24 May 2019
% Date modified: 7 March 2020 
%
%
function output = zeropad_odd_dimension(input,mode,matrixSize_o)

matrixSize = size(input);

% determine if a dimension needs to be zeropadded
padsize     = zeros(size(matrixSize));
for kd = 1:length(matrixSize)
    if mod(matrixSize(kd),2) == 1 && kd < 4
        padsize(kd) = 1;
    end
end

switch mode
    case 'pre'
        % zero padding if the dimension of the matrix is an odd number
        output = padarray(input, padsize, 0,'post');
        
    case 'post'
        % remove zero padding 
        output = input(1:matrixSize_o(1),1:matrixSize_o(2),1:matrixSize_o(3),:,:);
        
end

end