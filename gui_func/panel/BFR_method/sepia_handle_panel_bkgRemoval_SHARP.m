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
% Date last modified: 
%
%
function h = sepia_handle_panel_bkgRemoval_SHARP(hParent,h,position)

%% set default values
defaultRadius       = 4;
defaultthreshold    = 0.03;

%% Parent handle of SHARP panel children

h.bkgRemoval.panel.SHARP = uipanel(hParent,...
    'Title','SHARP',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','off');
    
%% children of SHARP panel

    % text|edit field pair: radius
    h.bkgRemoval.SHARP.text.radius = uicontrol('Parent',h.bkgRemoval.panel.SHARP,'Style','text',...
        'String','Radius (voxel):',...
        'units','normalized','position',[0.01 0.75 0.2 0.2],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Radius of spherical mean kernel');
    h.bkgRemoval.SHARP.edit.radius = uicontrol('Parent',h.bkgRemoval.panel.SHARP,'Style','edit',...
        'String',num2str(defaultRadius),...
        'units','normalized','position',[0.25 0.75 0.2 0.2],...0
        'backgroundcolor','white');

    % text|edit field pair: threshold
    h.bkgRemoval.SHARP.text.threshold = uicontrol('Parent',h.bkgRemoval.panel.SHARP,'Style','text',...
        'String','Threshold:',...
        'units','normalized','position',[0.01 0.5 0.2 0.2],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Threshold used in truncated SVD');
    h.bkgRemoval.SHARP.edit.threshold = uicontrol('Parent',h.bkgRemoval.panel.SHARP,'Style','edit',...
        'String',num2str(defaultthreshold),...
        'units','normalized','position',[0.25 0.5 0.2 0.2],...
        'backgroundcolor','white');

%% set callbacks
set(h.bkgRemoval.SHARP.edit.radius,   	'Callback', {@EditInputMinMax_Callback,defaultRadius,       1,0});
set(h.bkgRemoval.SHARP.edit.threshold,	'Callback', {@EditInputMinMax_Callback,defaultthreshold,    0,0});

end