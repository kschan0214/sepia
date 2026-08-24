%% bgColor = sepia_theme_bg_color()
%
% Output
% --------------
% bgColor   : background colour (RGB triplet) for interactive uicontrols
%             (edit/popup/pushbutton fields)
%
% Description: Returns a background colour for interactive uicontrols
% that keeps them legible under both MATLAB's classic light desktop and
% its dark desktop theme (dark theme support for figure/uicontrol-based
% apps was introduced in MATLAB R2025a).
%
% SEPIA never sets its GUI figure's 'Color' property explicitly, so on
% MATLAB releases with theme support the figure's resolved colour
% (get(gcf,'color')) already automatically tracks the active desktop
% theme. This function reuses that resolved colour as a robust,
% version-independent signal instead of querying any theme-specific API,
% which would not exist on pre-R2025a MATLAB. See also
% sepia_theme_fg_color.m for the matching text colour.
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 24 August 2026 (v1.3.0)
%
function bgColor = sepia_theme_bg_color()

try
    figColor = get(gcf,'color');
catch
    figColor = [0.94 0.94 0.94]; % classic MATLAB grey figure default
end

if mean(figColor) < 0.5
    % Dark desktop theme is active: use a shade a little lighter than the
    % figure/panel background so the control stays visually distinct.
    bgColor = min(figColor + 0.12, 1);
else
    % Classic light desktop theme (or a MATLAB release without dark
    % theme support): keep the original white background.
    bgColor = [1 1 1];
end

end
