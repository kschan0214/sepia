%% upsampled = fft_upsample_complex(img, newSize)
%
% Input
% --------------
% img           :   2D or 3D image
% newSize       :   target size
%
% Output
% --------------
% upsampled     : upsampled image
%
% Description: Upsample 2D or 3D data using k-space zeo-filling
%
% Kwok-shing Chan @ MGH
% kchan2@mgh.harvard.edu
% Date created: 4 July 2025
% Date modified: 
% 
function upsampled = fft_upsample_complex(img, newSize)

    % FFT to k-space
    kspace = fftshift(fftn(ifftshift(img)));

    % Get original size
    origSize = size(img);
    
    % Initialize zero-padded k-space
    kspace_padded = zeros(newSize, class(kspace));

    % Calculate start and end indices for placing original k-space
    startIdx    = floor((newSize - origSize)/2) + 1;
    endIdx      = startIdx + origSize - 1;

    % Fill in the center of the zero-padded array
    switch ndims(img)
        case 2
            kspace_padded(startIdx(1):endIdx(1), startIdx(2):endIdx(2)) = kspace;
        case 3
            kspace_padded(startIdx(1):endIdx(1), startIdx(2):endIdx(2), startIdx(3):endIdx(3)) = kspace;
        otherwise
            error('Only 2D or 3D images supported');
    end

    % Inverse FFT to image space
    upsampled = fftshift(ifftn(ifftshift(kspace_padded)));

    % Optional: scale to match original intensity
    upsampled = upsampled * prod(newSize) / prod(origSize);
end
