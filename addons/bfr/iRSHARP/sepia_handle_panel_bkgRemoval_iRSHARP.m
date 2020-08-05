%% h = sepia_handle_panel_qsm_iterTik(hParent,h,position)
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
% Description: This GUI function creates a panel for iterTik method
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 4 August 2018
% Date modified: 
%
%
function h = sepia_handle_panel_bkgRemoval_iRSHARP(hParent,h,position)

%% set default values
defaultRadius      	= 8;
defaultThreshold	= 0.05;
defaultConstant   	= 0.25;

%% Tooltips
tooltip.bkgRemoval.iRSHARP.radius   	= 'SHARP radius in mm';
tooltip.bkgRemoval.iRSHARP.threshold	= 'Truncated SVD threshold';
tooltip.bkgRemoval.iRSHARP.constant   	= 'Adjusting constant';

%% layout of the panel
nrow        = 4;
rspacing    = 0.03;
ncol        = 2;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent handle of iterTik panel children

h.bkgRemoval.panel.iRSHARP = uipanel(hParent,...
        'Title','iRSHARP',...
        'position',position,...
        'backgroundcolor',get(h.fig,'color'),'Visible','off');

%% Children of iterTik panel
    
    panelParent = h.bkgRemoval.panel.iRSHARP;

    % width of each element in a functional column, in normalised unit
    wratio = 0.5;
    
    % col 1, row 1
    % text|edit field pair: maximum iterations
    [h.bkgRemoval.iRSHARP.text.radius,h.bkgRemoval.iRSHARP.edit.radius] = sepia_construct_text_edit(...
        panelParent,'SHARP radius (mm):',      	defaultRadius,      [left(1) bottom(1) width height], wratio);
    
    % col 1, row 3
    % text|edit field pair: maximum iterations
    [h.bkgRemoval.iRSHARP.text.threshold,h.bkgRemoval.iRSHARP.edit.threshold] = sepia_construct_text_edit(...
        panelParent,'Truncated SVD threshold:',	defaultThreshold, 	[left(1) bottom(2) width height], wratio);
    
    % col 1, row 4
    % text|edit field pair: tolerance
    [h.bkgRemoval.iRSHARP.text.constant,h.bkgRemoval.iRSHARP.edit.constant] = sepia_construct_text_edit(...
        panelParent,'Constant:',                defaultConstant,   	[left(1) bottom(3) width height], wratio);
    
    
%% set tooltips
set(h.bkgRemoval.iRSHARP.text.radius,        	'Tooltip',tooltip.bkgRemoval.iRSHARP.radius);
set(h.bkgRemoval.iRSHARP.text.threshold,       	'Tooltip',tooltip.bkgRemoval.iRSHARP.threshold);
set(h.bkgRemoval.iRSHARP.text.constant,        	'Tooltip',tooltip.bkgRemoval.iRSHARP.constant);

%% set callbacks
set(h.bkgRemoval.iRSHARP.edit.radius,          'Callback', {@EditInputMinMax_Callback,defaultRadius,   0,0});
set(h.bkgRemoval.iRSHARP.edit.threshold,       'Callback', {@EditInputMinMax_Callback,defaultThreshold,0,0});
set(h.bkgRemoval.iRSHARP.edit.constant,        'Callback', {@EditInputMinMax_Callback,defaultConstant, 0,0});

end
