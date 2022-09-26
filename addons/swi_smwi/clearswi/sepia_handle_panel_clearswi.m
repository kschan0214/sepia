%% h = sepia_handle_panel_clearswi(hParent,h,position)
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
% Description: This GUI function creates a panel for CLEARSWI method
%
% Korbinian Eckstein @ UQ
% korbinian90@gmail.com
% Date created: 26 September 2022
%
%
function h = sepia_handle_panel_clearswi(hParent,h,position)

%% set default values
defaultPhaseScalingStrength = 4;
defaultFilterSize   = '[4,4,0]';
defaultEchoes = ':';
defaultEchoCombineMethodAdd = 1;
defaultmIP          = 4;

%% layout of the panel
nrow        = 4;
rspacing    = 0.01;
ncol        = 2;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent handle of panel children

h.swismwi.panel.clearswi = uipanel(hParent,...
        'Title','SWI',...
        'position',position,...
        'backgroundcolor',get(h.fig,'color'),'Visible','on');

%% Children of panel

    panelParent = h.swismwi.panel.clearswi;

    % width of each element in a functional column, in normalised unit
    wratio = 0.5;
   
    % row 1, col 1
    % text|edit field pair: Phase Function
    [h.swismwi.clearswi.text.phaseScalingType ,h.swismwi.clearswi.popup.phaseScalingType] = sepia_construct_text_popup(...
        panelParent,'Phase Function:', {'tanh', 'negativetanh', 'positive', 'negative', 'triangular'}, [left(1) bottom(1) width height], wratio);
   
    % row 2, col 1
    % text|edit field pair: Contrast
    [h.swismwi.clearswi.text.phaseScalingStrength ,h.swismwi.clearswi.edit.phaseScalingStrength] = sepia_construct_text_edit(...
        panelParent,'Phase Contrast:', defaultPhaseScalingStrength, [left(1) bottom(2) width height], wratio);

    
    % row 3, col 1
    % text|edit field pair: filter size
    [h.swismwi.clearswi.text.filterSize ,h.swismwi.clearswi.edit.filterSize] = sepia_construct_text_edit(...
        panelParent,'Filter size:', defaultFilterSize, [left(1) bottom(3) width height], wratio);
    
    % row 4, col 1
    % text|popup field pair: Unwrapping Algorithm
    [h.swismwi.clearswi.text.unwrappingAlgorithm ,h.swismwi.clearswi.popup.unwrappingAlgorithm] = sepia_construct_text_popup(...
        panelParent,'Unwrapping Algorithm:', {'laplacian', 'romeo', 'laplacianslice'}, [left(1) bottom(4) width height], wratio);
    
    % row 1, col 2
    % text|popup field pair: Magnitude Echo Combination
    [h.swismwi.clearswi.text.echoCombineMethod ,h.swismwi.clearswi.popup.echoCombineMethod] = sepia_construct_text_popup(...
        panelParent,'Magnitude Echo Combination:', {'SNR', 'average', 'echo', 'simulated echo'}, [left(2) bottom(1) width height], wratio);
    
    % row 2, col 2/1
    % text|edit field pair: Echoes to include
    [h.swismwi.clearswi.text.echoes ,h.swismwi.clearswi.edit.echoes] = sepia_construct_text_edit(...
        panelParent,'Echoes to include:', defaultEchoes, [left(2) bottom(2) width/2 height], wratio);
    
    % row 2, col 2/2
    % text|edit field pair: Additional Echo Combine Parameter
    [h.swismwi.clearswi.text.echoCombineMethodAdd ,h.swismwi.clearswi.edit.echoCombineMethodAdd] = sepia_construct_text_edit(...
        panelParent,'echonumber/echotime:', defaultEchoCombineMethodAdd, [left(2)+width/2 bottom(2) width/2 height], wratio);    
    
    % row 3, col 2/1
    % checkbox: Softplus Magnitude Scaling
    h.swismwi.clearswi.checkbox.softplusScaling = uicontrol('Parent',h.swismwi.panel.clearswi,'Style','checkbox',...
        'String','Softplus Magnitude Scaling',...
        'units','normalized','position',[left(2) bottom(3) width/2 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),'Value',true);
    % row 3, col 2/2
    % checkbox: Sensitivity Correction
    h.swismwi.clearswi.checkbox.sensitivityCorrection = uicontrol('Parent',h.swismwi.panel.clearswi,'Style','checkbox',...
        'String','Sensitivity Correction',...
        'units','normalized','position',[left(2)+width/2 bottom(3) width/2 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),'Value',true);
    % row 4, col 2
    % checkbox: save minimum intensity projection image
    [h.swismwi.clearswi.checkbox.mIP,h.swismwi.clearswi.edit.mIP] = sepia_construct_checkbox_edit(...
        panelParent,'Save mIP image, #slices', defaultmIP, [left(2) bottom(4) width height], wratio);
    set(h.swismwi.clearswi.checkbox.mIP,'Value',true);
    set(h.swismwi.clearswi.edit.mIP,    'enable','on');
    

%% set callbacks
set(h.swismwi.clearswi.checkbox.mIP,  	'Callback', {@CheckboxEditPair_Callback,h.swismwi.clearswi.edit.mIP,1});

set(h.swismwi.clearswi.edit.echoCombineMethodAdd, 'Enable', 'off'); % default: off
set(h.swismwi.clearswi.popup.echoCombineMethod,               'Callback', {@clearswi_echo_combine_selection_Callback,h});

%% Callback functions
function clearswi_echo_combine_selection_Callback(source,eventdata,h)
    if strcmp(source.String{source.Value,1}, 'echo') || strcmp(source.String{source.Value,1}, 'simulated echo')
        set(h.swismwi.clearswi.edit.echoCombineMethodAdd, 'Enable', 'on');
    else
        set(h.swismwi.clearswi.edit.echoCombineMethodAdd, 'Enable', 'off');
    end
end

end