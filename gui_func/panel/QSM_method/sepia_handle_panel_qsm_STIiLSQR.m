%% h = sepia_handle_panel_qsm_STIiLSQR(hParent,h,position)
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
% Description: This GUI function creates a panel for STIiLSQR method
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 1 June 2018
% Date modified: 3 April 2020 (v0.8.0)
%
%
function h = sepia_handle_panel_qsm_STIiLSQR(hParent,h,position)

%% set default values
defaultTol1         = 0.01;
defaultTol2         = 0.001;
defaultPadSize    	= 12;
defaultThreshold	= 0.01;
defaultMaxIter      = 100;

%% Tooltips
% tooltip.qsm.STIiLSQR.threshold	= '';
% tooltip.qsm.STIiLSQR.maxIter      = '';

%% layout of the panel
nrow        = 4;
rspacing    = 0.01;
ncol        = 2;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent handle of CFS panel children

h.qsm.panel.STIiLSQR = uipanel(hParent,...
    'Title','STI suite iLSQR',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','off');

%% Children of panel
    
    panelParent = h.qsm.panel.STIiLSQR;

    % width of each element in a functional column, in normalised unit
    wratio = 0.5;
    
    % col 1, row 1
    % text|edit field pair: threshold
    [h.qsm.STIiLSQR.text.threshold,h.qsm.STIiLSQR.edit.threshold] = sepia_construct_text_edit(...
        panelParent,'Threshold:', defaultThreshold, [left(1) bottom(1) width height], wratio);
   
    % col 1, row 2
	% text|edit field pair: maximum iterations
    [h.qsm.STIiLSQR.text.maxIter,h.qsm.STIiLSQR.edit.maxIter] = sepia_construct_text_edit(...
        panelParent,'Max. iterations:', defaultMaxIter, [left(1) bottom(2) width height], wratio);
    
    % col 1, row 3
	% text|edit field pair: tolerance 1
    [h.qsm.STIiLSQR.text.tol1,h.qsm.STIiLSQR.edit.tol1] = sepia_construct_text_edit(...
        panelParent,'Tolerance 1:', defaultTol1, [left(1) bottom(3) width height], wratio);
    
    % col 1, row 4
    % text|edit field pair: tolerance 2
    [h.qsm.STIiLSQR.text.tol2,h.qsm.STIiLSQR.edit.tol2] = sepia_construct_text_edit(...
        panelParent,'Tolerance 2:', defaultTol2, [left(1) bottom(4) width height], wratio);

    % col 2, row 1
    % text|edit field pair: pad size
    [h.qsm.STIiLSQR.text.padSize,h.qsm.STIiLSQR.edit.padSize] = sepia_construct_text_edit(...
        panelParent,'Pad size:', defaultPadSize, [left(2) bottom(1) width height], wratio);


%% set callbacks
set(h.qsm.STIiLSQR.edit.maxIter,	'Callback', {@EditInputMinMax_Callback,defaultMaxIter,  1,1});
set(h.qsm.STIiLSQR.edit.padSize, 	'Callback', {@EditInputMinMax_Callback,defaultPadSize,  1,0});
set(h.qsm.STIiLSQR.edit.threshold,	'Callback', {@EditInputMinMax_Callback,defaultThreshold,0,0});
set(h.qsm.STIiLSQR.edit.tol1,     	'Callback', {@EditInputMinMax_Callback,defaultTol1,     0,0});
set(h.qsm.STIiLSQR.edit.tol2,    	'Callback', {@EditInputMinMax_Callback,defaultTol2,     0,0});

end