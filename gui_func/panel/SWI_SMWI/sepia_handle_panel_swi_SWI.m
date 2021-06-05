%% h = sepia_handle_panel_qsm_FANSI(hParent,h,position)
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
% Description: This GUI function creates a panel for SWI method
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 1 June 2018
% Date modified: 4 April 2020 (v0.8.0)
%
%
function h = sepia_handle_panel_swi_SWI(hParent,h,position)

%% set default values
defaultM            = 4;
defaultthreshold    = 'pi';
defaultFilterSize   = 12;
defaultmIP          = 4;

%% layout of the panel
nrow        = 4;
rspacing    = 0.01;
ncol        = 2;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent handle of panel children

h.swi.panel.SWI = uipanel(hParent,...
        'Title','SWI',...
        'position',position,...
        'backgroundcolor',get(h.fig,'color'),'Visible','on');

%% Children of panel

    panelParent = h.swi.panel.SWI;

    % width of each element in a functional column, in normalised unit
    wratio = 0.5;
    
    % row 1, col 1
    % text|edit field pair: radius
    [h.swi.SWI.text.m ,h.swi.SWI.edit.m] = sepia_construct_text_edit(...
        panelParent,'Contrast:', defaultM, [left(1) bottom(1) width height], wratio);

    % row 2, col 1
    % text|edit field pair: threshold
    [h.swi.SWI.text.threshold ,h.swi.SWI.edit.threshold] = sepia_construct_text_edit(...
        panelParent,'Threshold (rad):', defaultthreshold, [left(1) bottom(2) width height], wratio);
   
    % row 3, col 1
    % text|edit field pair: filter size
    [h.swi.SWI.text.filterSize ,h.swi.SWI.edit.filterSize] = sepia_construct_text_edit(...
        panelParent,'Filter size:', defaultFilterSize, [left(1) bottom(3) width height], wratio);
    
    % row 4, col 1
    % text|popup field pair: method
    [h.swi.SWI.text.method ,h.swi.SWI.popup.method] = sepia_construct_text_popup(...
        panelParent,'Method:', {'default','multi-echo (testing)'}, [left(1) bottom(4) width height], wratio);
    
    % row 1, col 2
    % checkbox: save positive phase image
    h.swi.SWI.checkbox.positive = uicontrol('Parent',h.swi.panel.SWI,'Style','checkbox',...
        'String','Save positive phase weighted images',...
        'units','normalized','position',[left(2) bottom(1) width height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),'Value',true');
    % row 2, col 2
    % checkbox: save nagative phase image
    h.swi.SWI.checkbox.negative = uicontrol('Parent',h.swi.panel.SWI,'Style','checkbox',...
        'String','Save negative phase weighted images',...
        'units','normalized','position',[left(2) bottom(2) width height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'));
    % row 3, col 2
    % checkbox: save minimum intensity projection image
    [h.swi.SWI.checkbox.mIP,h.swi.SWI.edit.mIP] = sepia_construct_checkbox_edit(...
        panelParent,'Save mIP image, #slices', defaultmIP, [left(2) bottom(3) width height], wratio);
    set(h.swi.SWI.checkbox.mIP,'Value',true);
    set(h.swi.SWI.edit.mIP,    'enable','on');
    

%% set callbacks
set(h.swi.SWI.edit.m,           'Callback', {@EditInputMinMax_Callback,defaultM,   1,0});
set(h.swi.SWI.edit.threshold,	'Callback', {@EditInputMinMax_Callback,str2double(defaultthreshold),  0,0});
set(h.swi.SWI.edit.filterSize,  'Callback', {@EditInputMinMax_Callback,defaultFilterSize,       1,1});
set(h.swi.SWI.checkbox.mIP,     'Callback', {@CheckboxEditPair_Callback,h.swi.SWI.edit.mIP,1});
end