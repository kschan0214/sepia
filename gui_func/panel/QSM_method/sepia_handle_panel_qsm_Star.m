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
    
    % width of each element in a functional column, in normalised unit of
    % the functional column width
    subwidth(1) = width*0.5;
    subwidth(2) = width-subwidth(1);
    
    % text|edit field pair: pad size
    h.qsm.Star.text.padSize = uicontrol('Parent',h.qsm.panel.Star,...
        'Style','text',...
        'String','Pad size:',...
        'units','normalized','position',[left(1) bottom(1) subwidth(1) height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'Tooltip',tooltip.qsm.Star.padSize);
    h.qsm.Star.edit.padSize = uicontrol('Parent',h.qsm.panel.Star,...
        'Style','edit',...
        'String',num2str(defaultPadSize),...
        'units','normalized','position',[left(1)+subwidth(1) bottom(1) subwidth(2) height],...
        'backgroundcolor','white');

%% set callbacks
set(h.qsm.Star.edit.padSize, 'Callback', {@EditInputMinMax_Callback,defaultPadSize,1,0});

end