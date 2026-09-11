%% icon = sepia_theme_open_icon()
%
% Output
% --------------
% icon      : 16x16 CData image matrix for the folder/file 'open' icon
%             used on the GUI's push buttons
%
% Description: Returns the folder icon used for file/directory 'open'
% push buttons, automatically switching between the light-theme
% ('folder@0,3x.jpg') and dark-theme ('folder_dark@0,3x.jpg') variant
% depending on whether MATLAB's dark desktop theme is active (see
% sepia_theme_bg_color.m for how that is detected).
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 24 August 2026 (v1.3.0)
%
function icon = sepia_theme_open_icon()

try
    figColor = get(gcf,'color');
catch
    figColor = [0.94 0.94 0.94]; % classic MATLAB grey figure default
end

if mean(figColor) < 0.5
    iconFilename = 'folder_dark@0,3x.jpg';
else
    iconFilename = 'folder@0,3x.jpg';
end

icon = imread(iconFilename);
icon = imresize(icon,[1 1]*20);

end
