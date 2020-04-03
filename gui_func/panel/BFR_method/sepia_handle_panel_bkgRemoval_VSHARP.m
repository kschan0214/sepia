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
%
%
function h = sepia_handle_panel_bkgRemoval_VSHARP(hParent,h,position)

%% set default values
defaultMaxRadius = 10;
defaultMinRadius = 3;

%% Tooltips
tooltip.bkgRemoval.VSHARP.maxRadius = 'Maximum radius of spherical mean value kernel';
tooltip.bkgRemoval.VSHARP.minRadius = 'Minimum radius of spherical mean value kernel';

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

    % width of each element in a functional column, in normalised unit of
    % the functional column width
    subwidth(1) = width*0.5;
    subwidth(2) = width-subwidth(1);
    
    % row 1
    % text|edit field pair: maximum radius
    h.bkgRemoval.VSHARP.text.maxRadius = uicontrol('Parent',h.bkgRemoval.panel.VSHARP,...
        'Style','text',...
        'String','Max. radius (voxel):',...
        'units','normalized','position',[left(1) bottom(1) subwidth(1) height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip',tooltip.bkgRemoval.VSHARP.maxRadius);
    h.bkgRemoval.VSHARP.edit.maxRadius = uicontrol('Parent',h.bkgRemoval.panel.VSHARP,...
        'Style','edit',...
        'String',num2str(defaultMaxRadius),...
        'units','normalized','position',[left(1)+subwidth(1) bottom(1) subwidth(2) height],...
        'backgroundcolor','white');
    
    % row 2
    % text|edit field pair: minimum radius
    h.bkgRemoval.VSHARP.text.minRadius = uicontrol('Parent',h.bkgRemoval.panel.VSHARP,...
        'Style','text',...
        'String','Min. radius (voxel):',...
        'units','normalized','position',[left(1) bottom(2) subwidth(1) height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip',tooltip.bkgRemoval.VSHARP.minRadius);
    h.bkgRemoval.VSHARP.edit.minRadius = uicontrol('Parent',h.bkgRemoval.panel.VSHARP,...
        'Style','edit',...
        'String',num2str(defaultMinRadius),...
        'units','normalized','position',[left(1)+subwidth(1) bottom(2) subwidth(2) height],...
        'backgroundcolor','white');

%% set callbacks
set(h.bkgRemoval.VSHARP.edit.minRadius, 'Callback', {@EditVSHARPRadius_Callback,h});
set(h.bkgRemoval.VSHARP.edit.maxRadius, 'Callback', {@EditVSHARPRadius_Callback,h});

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
