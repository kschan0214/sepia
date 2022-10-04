%% h = sepia_handle_panel_qsm_CFS(hParent,h,position)
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
% Description: This GUI function creates a panel for CFS method
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 1 June 2018
% Date modified: 4 April 2020
%
%
function h = sepia_handle_panel_R2s_NLLS(hParent,h,position)

%% set default values
menuFitType = {'Magnitude','Complex','Mixed'};

%% Tooltips

%% layout of the panel
nrow        = 4;
rspacing    = 0.01;
ncol        = 2;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent handle of CFS panel children

h.r2s.panel.nlls = uipanel(hParent,...
    'Title','Trapezoidal',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','off');
        
%% Children of CFS panel
    
    panelParent = h.r2s.panel.nlls;

    % width of each element in a functional column, in normalised unit
    wratio = 0.5;
    
    % row 1, col 1
    % text|edit field pair: S0 parameter
    [h.r2s.nlls.text.fit,h.r2s.nlls.popup.fit] = sepia_construct_text_popup(...
        panelParent,'Fitting type:', menuFitType, [left(1) bottom(1) width height], wratio);

    h.r2s.nlls.checkbox.isParallel = uicontrol('Parent',panelParent,...
        'Style','checkbox',...
        'String','Enable parallel computing (parfor)',...
        'units','normalized','position',[left(1) bottom(2) width height]);

%% set tooltips
% set(h.r2s.trapezoidal.text.s0,      'Tooltip',tooltip.r2s.trapezoidal.s0);

%% set callbacks


end