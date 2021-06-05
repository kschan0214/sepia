%% h = sepia_handle_panel_bkgRemoval_SHARP(hParent,h,position)
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
% Date modified: 4 April 2020 (v0.8.0)
%
%
function h = sepia_handle_panel_swi_SMWI(hParent,h,position)

%% set default values
defaultM            = 4;
defaultthreshold    = 1;
defaultmIP          = 4;

%% layout of the panel
nrow        = 4;
rspacing    = 0.01;
ncol        = 2;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent handle of SMWI panel children

h.swi.panel.SMWI = uipanel(hParent,...
    'Title','SMWI',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','off');
    
%% children of SMWI panel

    panelParent = h.swi.panel.SMWI;

    % width of each element in a functional column, in normalised unit
    wratio = 0.5;

    % row 1, col 1
    % text|edit field pair: radius
    [h.swi.SMWI.text.m,h.swi.SMWI.edit.m] = sepia_construct_text_edit(...
        panelParent,'Contrast:', defaultM, [left(1) bottom(1) width height], wratio);

    % row 2, col 1
    % text|edit field pair: threshold
    [h.swi.SMWI.text.threshold,h.swi.SMWI.edit.threshold] = sepia_construct_text_edit(...
        panelParent,'Threshold (ppm):', defaultthreshold, [left(1) bottom(2) width height], wratio);
    
    % row 1, col 2
    % checkbox: save paramagetic image
    h.swi.SMWI.checkbox.paramagnetic = uicontrol('Parent',h.swi.panel.SMWI,'Style','checkbox',...
        'String','Save paramagnetic weighted images',...
        'units','normalized','position',[left(2) bottom(1) width height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),'Value',true');
    
    % row 2, col 2
    % checkbox: save diamagetic image
    h.swi.SMWI.checkbox.diamagnetic = uicontrol('Parent',h.swi.panel.SMWI,'Style','checkbox',...
        'String','Save diamagnetic weighted images',...
        'units','normalized','position',[left(2) bottom(2) width height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'));
    
    % row 3, col 2
    % checkbox: save nagative phase image
    [h.swi.SMWI.checkbox.mIP,h.swi.SMWI.edit.mIP] = sepia_construct_checkbox_edit(...
        panelParent,'Save mIP image, #slices', defaultmIP, [left(2) bottom(3) width height], wratio);
    set(h.swi.SMWI.checkbox.mIP,'Value',true);
    set(h.swi.SMWI.edit.mIP,    'enable','on');

%% Tooltip
% set(h.swi.SMWI.text.m,          'Tooltip', '')
% set(h.swi.SMWI.text.threshold,  'Tooltip', '')

%% set callbacks
set(h.swi.SMWI.edit.m,          'Callback', {@EditInputMinMax_Callback,defaultM,     	1,0});
set(h.swi.SMWI.edit.threshold,	'Callback', {@EditInputMinMax_Callback,defaultthreshold,0,0});
set(h.swi.SMWI.checkbox.mIP,    'Callback', {@CheckboxEditPair_Callback,h.swi.SMWI.edit.mIP,1});
end