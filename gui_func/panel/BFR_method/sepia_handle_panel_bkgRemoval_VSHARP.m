%% h = sepia_handle_panel_bkgRemoval_VSHARP(hParent,h,position)
%
% Input
% --------------
% hParent       : parent handle of this panel
% h             : global structure contains all handles
% position      : position of this panel
%
% Output
% --------------
% h             : global structure contains all new and other handles
%
% Description: This GUI function creates a panel for VSHARP 
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 1 June 2018
% Date modified: 4 April 2020 (v0.8.0)
% Date modified: 23 August 2026 (v1.3.0)
%
function h = sepia_handle_panel_bkgRemoval_VSHARP(hParent,h,position)

%% set default values
defaultMaxRadius = 12;  % mm
defaultMinRadius = 1;   % mm
defaultThreshold = 0.05;

%% Tooltips
tooltip.bkgRemoval.VSHARP.maxRadius = 'Maximum radius of spherical mean value kernel, in mm';
tooltip.bkgRemoval.VSHARP.minRadius = 'Minimum radius of spherical mean value kernel, in mm';
tooltip.bkgRemoval.VSHARP.threshold = 'Threshold used in the k-space deconvolution (Truncated SVD)';

%% layout of the panel
nrow        = 4;
rspacing    = 0.01;
ncol        = 2;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent handle of VSHARP panel children

h.bkgRemoval.panel.VSHARP = uipanel(hParent,...
    'Title','Variable SHARP',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','off');

%% Children of VSHARP panel

    panelParent = h.bkgRemoval.panel.VSHARP;

    % width of each element in a functional column, in normalised unit
    wratio = 0.5;
    
    % row 1, col 1
    % text|edit field pair: maximum radius
    [h.bkgRemoval.VSHARP.text.maxRadius,h.bkgRemoval.VSHARP.edit.maxRadius] = sepia_construct_text_edit(...
        panelParent,'Max. radius (mm):', defaultMaxRadius, [left(1) bottom(1) width height], wratio);

    % row 2, col 1
    % text|edit field pair: minimum radius
    [h.bkgRemoval.VSHARP.text.minRadius,h.bkgRemoval.VSHARP.edit.minRadius] = sepia_construct_text_edit(...
        panelParent,'Min. radius (mm):', defaultMinRadius, [left(1) bottom(2) width height], wratio);

    % row 3, col 1
    % text|edit field pair: threshold
    [h.bkgRemoval.VSHARP.text.threshold,h.bkgRemoval.VSHARP.edit.threshold] = sepia_construct_text_edit(...
        panelParent,'Threshold:', defaultThreshold, [left(1) bottom(3) width height], wratio);


%% set tooltips
set(h.bkgRemoval.VSHARP.text.maxRadius,	'Tooltip',tooltip.bkgRemoval.VSHARP.maxRadius);
set(h.bkgRemoval.VSHARP.text.minRadius,	'Tooltip',tooltip.bkgRemoval.VSHARP.minRadius);
set(h.bkgRemoval.VSHARP.text.threshold,	'Tooltip',tooltip.bkgRemoval.VSHARP.threshold);

%% set callbacks
set(h.bkgRemoval.VSHARP.edit.minRadius, 'Callback', {@EditVSHARPRadius_Callback,h});
set(h.bkgRemoval.VSHARP.edit.maxRadius, 'Callback', {@EditVSHARPRadius_Callback,h});
set(h.bkgRemoval.VSHARP.edit.threshold, 'Callback', {@EditInputMinMax_Callback,defaultThreshold,0,0});

end

%% Callback
function EditVSHARPRadius_Callback(source,eventdata,h)
% constraint the minimum of maximum radius is always larger then the
% minimum radius

% global h

% check minimum of minimum radius input
if str2double(h.bkgRemoval.VSHARP.edit.minRadius.String)<0
    h.bkgRemoval.VSHARP.edit.minRadius.String = num2str(0);
end

% if the minimum radius is not integer then rounds it to interger
h.bkgRemoval.VSHARP.edit.minRadius.String = num2str(round(str2double(h.bkgRemoval.VSHARP.edit.minRadius.String)));

% ensure maximum radius is always larger then minimum radius
if str2double(h.bkgRemoval.VSHARP.edit.maxRadius.String) <= str2double(h.bkgRemoval.VSHARP.edit.minRadius.String)
    h.bkgRemoval.VSHARP.edit.maxRadius.String = num2str(str2double(h.bkgRemoval.VSHARP.edit.minRadius.String) +1);
end

% if the maximum radius is not integer then rounds it to interger
h.bkgRemoval.VSHARP.edit.maxRadius.String = num2str(round(str2double(h.bkgRemoval.VSHARP.edit.maxRadius.String)));

end
