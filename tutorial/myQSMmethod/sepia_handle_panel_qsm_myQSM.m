%% h = sepia_handle_panel_qsm_TKD(hParent,h,position)
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
function h = sepia_handle_panel_qsm_myQSM(hParent,h,position)

%% set default values
defaultThreshold = 0.15;

%% Tooltips
tooltip.qsm.myQSM.threshold	= 'K-space threshold';

%% layout of the panel
nrow        = 4;
rspacing    = 0.01;
ncol        = 2;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent handle of TKD panel children

h.qsm.panel.myQSM = uipanel(hParent,...
    'Title','Thresholded k-space division (TKD)',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','on');

%% Children of TKD panel
    
    panelParent = h.qsm.panel.myQSM;

    % width of each element in a functional column, in normalised unit
    wratio = 0.5;
    
    % row 1, col 1
    % text|edit field pair: threshold
    [h.qsm.myQSM.text.threshold,h.qsm.myQSM.edit.threshold] = sepia_construct_text_edit(...
        panelParent,'Threshold (0-1):', defaultThreshold, [left(1) bottom(1) width height], wratio);


%% set tooltips
set(h.qsm.myQSM.text.threshold, 'Tooltip',tooltip.qsm.myQSM.threshold);

%% set callbacks
set(h.qsm.myQSM.edit.threshold, 'Callback', {@EditInputMinMax_Callback,defaultThreshold,0,0,1});

end