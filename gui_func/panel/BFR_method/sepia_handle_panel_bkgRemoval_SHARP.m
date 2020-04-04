%% h = sepia_handle_panel_bkgRemoval_SHARP(hParent,h,position)
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
% Description: This GUI function creates a panel for LBV method
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 1 June 2018
% Date modified: 3 April 2020
%
%
function h = sepia_handle_panel_bkgRemoval_SHARP(hParent,h,position)

%% set default values
defaultRadius       = 4;
defaultthreshold    = 0.03;

%% Tooltips
tooltip.bkgRemoval.SHARP.radius     = 'Radius of spherical mean value kernel';
tooltip.bkgRemoval.SHARP.threshold	= 'Threshold used in Truncated SVD';

%% layout of the panel
nrow        = 4;
rspacing    = 0.01;
ncol        = 2;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent handle of SHARP panel children

h.bkgRemoval.panel.SHARP = uipanel(hParent,...
    'Title','SHARP',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','off');
    
%% children of SHARP panel

    panelParent = h.bkgRemoval.panel.SHARP;

    % width of each element in a functional column, in normalised unit
    wratio = 0.5;
    
    % row 1, col 1
    % text|edit field pair: radius
    [h.bkgRemoval.SHARP.text.radius,h.bkgRemoval.SHARP.edit.radius] = sepia_construct_text_edit(...
        panelParent,'SMV radius (voxel):',	defaultRadius,      [left(1) bottom(1) width height], wratio);

    % row 2
    % text|edit field pair: threshold
    [h.bkgRemoval.SHARP.text.threshold,h.bkgRemoval.SHARP.edit.threshold] = sepia_construct_text_edit(...
        panelParent,'Threshold:',         	defaultthreshold,	[left(1) bottom(2) width height], wratio);
    

%% set tooltips
set(h.bkgRemoval.SHARP.text.radius,     'Tooltip',tooltip.bkgRemoval.SHARP.radius);
set(h.bkgRemoval.SHARP.text.threshold,	'Tooltip',tooltip.bkgRemoval.SHARP.threshold);

%% set callbacks
set(h.bkgRemoval.SHARP.edit.radius,   	'Callback', {@EditInputMinMax_Callback,defaultRadius,       1,0});
set(h.bkgRemoval.SHARP.edit.threshold,	'Callback', {@EditInputMinMax_Callback,defaultthreshold,    0,0});

end