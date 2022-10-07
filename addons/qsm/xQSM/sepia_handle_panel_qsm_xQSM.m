%% h = sepia_handle_panel_qsm_xQSM(hParent,h,position)
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
function h = sepia_handle_panel_qsm_xQSM(hParent,h,position)

%% set default values
menuSolver       = {'xQSM_invivo', 'xQSM_syn', 'xQSM_invivo_withNoiseLayer', 'Unet_invivo', 'Unet_syn'};

%% Tooltips
tooltip.qsm.xQSM.solver   	 = 'Select a Pre-trained network';

%% layout of the panel
nrow        = 4;
rspacing    = 0.01;
ncol        = 2;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent handle of CFS panel children

h.qsm.panel.xQSM = uipanel(hParent,...
    'Title','xQSM',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','off');

%% Children of CFS panel

    panelParent = h.qsm.panel.xQSM;

    % width of each element in a functional column, in normalised unit
    wratio = 0.5;
    
    % col 1, row 1
    % text|edit field pair: tolerance
    [h.qsm.xQSM.text.solver,h.qsm.xQSM.popup.solver] = sepia_construct_text_popup(...
        panelParent,'Solver:',                  menuSolver,         [left(1) bottom(1) width height], wratio);

    
%% set tooltips
set(h.qsm.xQSM.text.solver,         'Tooltip',tooltip.qsm.xQSM.solver);

%% set callbacks

end