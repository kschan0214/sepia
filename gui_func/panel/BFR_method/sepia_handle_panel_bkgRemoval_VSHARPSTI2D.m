%% h = sepia_handle_panel_bkgRemoval_VSHARPSTI(hParent,h,position)
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
% Description: This GUI function creates a panel for STI suite VSHARP 
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 1 June 2018
% Date modified: 3 April 2020 (v0.8.0)
%
%
function h = sepia_handle_panel_bkgRemoval_VSHARPSTI2D(hParent,h,position)

defaultSMVSize = 12;

%% Tooltips
tooltip.bkgRemoval.VSHARPSTI2D.smvSize = 'Spherical mean value kernel size';

%% layout of the panel
nrow        = 4;
rspacing    = 0.01;
ncol        = 2;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent handle of VSHARPSTI panel children

h.bkgRemoval.panel.VSHARPSTI2D = uipanel(hParent,...
    'Title','STI suite Variable SHARP (2D) for 2D EPI corrected background field removal',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','off');

%% Child of VSHARPSTI panel

    panelParent = h.bkgRemoval.panel.VSHARPSTI2D;

    % width of each element in a functional column, in normalised unit
    wratio = 0.5;
    
    % row 1, col 1
    % text|edit field pair: SMV size
    [h.bkgRemoval.VSHARPSTI2D.text.smvSize,h.bkgRemoval.VSHARPSTI2D.edit.smvSize] = sepia_construct_text_edit(...
        panelParent,'SMV size (mm):', defaultSMVSize, [left(1) bottom(1) width height], wratio);
    

%% set tooltips
set(h.bkgRemoval.VSHARPSTI2D.text.smvSize,	'Tooltip', tooltip.bkgRemoval.VSHARPSTI2D.smvSize);

%% set callbacks
set(h.bkgRemoval.VSHARPSTI2D.edit.smvSize,    'Callback', {@EditInputMinMax_Callback,defaultSMVSize,1,0});

end