%% h = sepia_handle_panel_qsm_Star(hParent,h,position)
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
% Description: This GUI function creates a panel for FANSI method
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 1 June 2018
% Date modified: 3 April 2020 (v0.8.0)
%
%
function h = sepia_handle_panel_qsm_Star(hParent,h,position)

%% set default values
defaultPadSize      = 12;

%% Tooltips
tooltip.qsm.Star.padSize   	= 'Zero padding size';

%% layout of the panel
nrow        = 4;
rspacing    = 0.01;
ncol        = 2;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent handle of CFS panel children

h.qsm.panel.Star = uipanel(hParent,...
    'Title','Star-QSM',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','off');

%% Children of CFS panel

    panelParent = h.qsm.panel.Star;

    % width of each element in a functional column, in normalised unit
    wratio = 0.5;

    % col 1, row 1
    % text|edit field pair: pad size
    [h.qsm.Star.text.padSize,h.qsm.Star.edit.padSize] = sepia_construct_text_edit(...
        panelParent,'Pad size:', defaultPadSize, [left(1) bottom(1) width height], wratio);
    

%% set tooltips
set(h.qsm.Star.text.padSize, 'Tooltip', tooltip.qsm.Star.padSize);

%% set callbacks
set(h.qsm.Star.edit.padSize, 'Callback', {@EditInputMinMax_Callback,defaultPadSize,1,0});

end