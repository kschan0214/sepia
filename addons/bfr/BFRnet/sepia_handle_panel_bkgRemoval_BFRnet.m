%% h = sepia_handle_panel_bkgRemoval_BFRnet(hParent,h,position)
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
% Description: This GUI function creates a panel for NDI method
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 24 September 2022
% Date modified:
%
%
function h = sepia_handle_panel_bkgRemoval_BFRnet(hParent,h,position)

%% set default values

%% Tooltips

%% layout of the panel
nrow        = 4;
rspacing    = 0.01;
ncol        = 2;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent handle of BFRnet panel children

h.bkgRemoval.panel.BFRnet = uipanel(hParent,...
    'Title','BFRnet',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','off');

%% Children of BFRnet panel

%% set tooltips

end