%% fgColor = sepia_theme_fg_color()
%
% Output
% --------------
% fgColor   : text (foreground) colour (RGB triplet) for interactive
%             uicontrols (edit/popup/pushbutton fields)
%
% Description: Companion function to sepia_theme_bg_color.m - returns the
% text colour that keeps typed/selected text legible against the colour
% returned by sepia_theme_bg_color.m, under both MATLAB's classic light
% desktop and its dark desktop theme (introduced in MATLAB R2025a).
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 24 August 2026 (v1.3.0)
%
function fgColor = sepia_theme_fg_color()

try
    figColor = get(gcf,'color');
catch
    figColor = [0.94 0.94 0.94]; % classic MATLAB grey figure default
end

if mean(figColor) < 0.5
    % Dark desktop theme is active: white text on the lightened-dark
    % background returned by sepia_theme_bg_color.m
    fgColor = [1 1 1];
else
    % Classic light desktop theme: black text on white
    fgColor = [0 0 0];
end

end
