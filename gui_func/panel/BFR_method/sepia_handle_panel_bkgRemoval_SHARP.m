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

    % width of each element in a functional column, in normalised unit of
    % the functional column width
    subwidth(1) = width*0.5;
    subwidth(2) = width-subwidth(1);
    
    % row 1
    % text|edit field pair: radius
    h.bkgRemoval.SHARP.text.radius = uicontrol('Parent',h.bkgRemoval.panel.SHARP,'Style','text',...
        'String','Radius (voxel):',...
        'units','normalized','position',[left(1) bottom(1) subwidth(1) height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip',tooltip.bkgRemoval.SHARP.radius);
    h.bkgRemoval.SHARP.edit.radius = uicontrol('Parent',h.bkgRemoval.panel.SHARP,'Style','edit',...
        'String',num2str(defaultRadius),...
        'units','normalized','position',[left(1)+subwidth(1) bottom(1) subwidth(2) height],...
        'backgroundcolor','white');

    % row 2
    % text|edit field pair: threshold
    h.bkgRemoval.SHARP.text.threshold = uicontrol('Parent',h.bkgRemoval.panel.SHARP,'Style','text',...
        'String','Threshold:',...
        'units','normalized','position',[left(1) bottom(2) subwidth(1) height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Threshold used in truncated SVD');
    h.bkgRemoval.SHARP.edit.threshold = uicontrol('Parent',h.bkgRemoval.panel.SHARP,'Style','edit',...
        'String',num2str(defaultthreshold),...
        'units','normalized','position',[left(1)+subwidth(1) bottom(2) subwidth(2) height],...
        'backgroundcolor','white',...
        'Tooltip',tooltip.bkgRemoval.SHARP.threshold);

%% set callbacks
set(h.bkgRemoval.SHARP.edit.radius,   	'Callback', {@EditInputMinMax_Callback,defaultRadius,       1,0});
set(h.bkgRemoval.SHARP.edit.threshold,	'Callback', {@EditInputMinMax_Callback,defaultthreshold,    0,0});

end