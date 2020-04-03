%% h = sepia_handle_panel_phaseUnwrap(hParent,h,position)
%
% Input
% --------------
% hParent       : parent handle of this panel
% hFig          : handle of the GUI
% h             : global structure contains all handles
% position      : starting position of this panel
%
% Output
% --------------
% h             : global structure contains all new and other handles
%
% Description: This GUI function creates a panel for phase unwrapping method control
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 16 April 2018
% Date modified: 18 April 2018
% Date modified: 12 June 2018
% Date modified: 3 April 2020 (v0.8.0)
%
%
function h = sepia_handle_panel_phaseUnwrap(hParent,h,position)

% load universal variable for methods' names
sepia_universal_variables;

% set default value
defaultThreshold = 0.5;

% tooltips
tooltip.unwrap.panel.method     = 'Select a method to combine field maps from multi-echo data';
tooltip.unwrap.panel.unwrap     = 'Select a phase unwrapping algorithm for spatial unwrapping';
tooltip.unwrap.panel.eddy       = 'Correct phase offset between odd and even number echoes if bipolar readout was used in acquisition';
tooltip.unwrap.panel.exclude    = ['Apply threshold on relative residual to exclude unreliable voxels based on ',...
                                    'mono-exponential decay model (only available for non-Laplacian methods)'];
tooltip.unwrap.panel.saveUnwrap = ['If total field is computed using optimal weights combination with multi-echo data, ',...
                                   'it is possible to save the unwrapped phase for each echo.'];

%% layout of the panel
nrow        = 5;
rspacing    = 0.02;
ncol        = 2;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% design of this panel
% Parent handle of phase unwrapping panel
h.StepsPanel.phaseUnwrap = uipanel(hParent,...
    'Title','Total field recovery and phase unwrapping',...
    'fontweight', 'bold',...
    'position',[position(1) position(2) 0.95 0.2]);
    
    % width of each element in a functional column, in normalised unit of
    % the functional column width
    subwidth(1) = 0.5;
    subwidth(2) = 1-subwidth(1);
    % Temporo-spatial unwrapping methods, 'text|popup' functional column,
    % 1st row
    h.phaseUnwrap.text.phaseCombMethod = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,...
        'Style','text','String','Echo phase combination:',...
        'units','normalized','position',[left(1) bottom(1) width*subwidth(1) height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip',tooltip.unwrap.panel.method);
    h.phaseUnwrap.popup.phaseCombMethod = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,...
        'Style','popup',...
        'String',methodEchoCombineName,...
        'units','normalized','position',[left(1)+width*subwidth(1) bottom(1) width*subwidth(2) height]); 
    
    % phase unwrapping method, 'text|popup' functional column, 2nd row
    h.phaseUnwrap.text.phaseUnwrap = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,...
        'Style','text','String','Phase unwrapping:',...
        'units','normalized','position',[left(1) bottom(2) width*subwidth(1) height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip',tooltip.unwrap.panel.unwrap);
    h.phaseUnwrap.popup.phaseUnwrap = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,...
        'Style','popup',...
        'String',methodUnwrapName,...
        'units','normalized','position',[left(1)+width*subwidth(1) bottom(2) width*subwidth(2) height]); 

    % eddy current correction for bipolar readout, 'checkbox' functional
    % column, 3rd row
    h.phaseUnwrap.checkbox.eddyCorrect = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,...
        'Style','checkbox','String','Bipolar readout correction',...
        'units','normalized','Position',[0.01 bottom(3) 1 height],...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip',tooltip.unwrap.panel.eddy);
    
    % exclusion of unreliable voxels, 'checkbox|field|text|popup', 4th row
    h.phaseUnwrap.checkbox.excludeMask = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,...
        'Style','checkbox','String','Exclude voxels using residual, threshold:',...
        'units','normalized','position',[0.01 bottom(4) 0.3 height],...
        'backgroundcolor',get(h.fig,'color'),...
        'Enable','off',...
        'tooltip',tooltip.unwrap.panel.exclude);
    h.phaseUnwrap.edit.excludeMask = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,...
        'Style','edit',...
        'String',num2str(defaultThreshold),...
        'units','normalized','position',[0.31 bottom(4) 0.04 height],...
        'backgroundcolor','white',...
        'Enable','off',...
        'tooltip','Higher value means accepting larger error between the data fitting and measurement');
    % excluding method
    h.phaseUnwrap.text.excludeMethod = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,...
        'Style','text','String','and apply in ',...
        'units','normalized','position',[0.36 bottom(4) 0.1 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'));
    h.phaseUnwrap.popup.excludeMethod = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,...
        'Style','popup',...
        'String',methodExcludedName,...
        'Enable','off',...
        'units','normalized','position',[0.46 bottom(4) 0.23 height]); 
    
    % save unwrapped echo phase option, 'checkbox', 5th row
    h.phaseUnwrap.checkbox.saveEchoPhase = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,...
        'Style','checkbox','String','Save unwrapped echo phase',...
        'units','normalized','position',[0.01 bottom(5) 0.3 height],...
        'backgroundcolor',get(h.fig,'color'),...
        'Enable','on',...
        'tooltip',tooltip.unwrap.panel.saveUnwrap);

%% set callback functions
set(h.phaseUnwrap.checkbox.excludeMask,	'Callback', {@CheckboxEditPair_Callback,{h.phaseUnwrap.edit.excludeMask,h.phaseUnwrap.popup.excludeMethod},1});
set(h.phaseUnwrap.edit.excludeMask, 	'Callback', {@EditInputMinMax_Callback,defaultThreshold,0,0,1});
set(h.phaseUnwrap.popup.phaseUnwrap,    'Callback', {@popupPhaseUnwrap_Callback,h});
set(h.phaseUnwrap.popup.phaseCombMethod,'Callback', {@popupPhaseCombine_Callback,h});
        
end

%% Callback functions
% Phase unwrapping method specific panel setting
function popupPhaseUnwrap_Callback(source,eventdata,h)

sepia_universal_variables;

% get selected background removal method
method = source.String{source.Value,1} ;

% method the user chosen will affect if exclusion method can be used or not 
switch method
    case methodUnwrapName{1}    % Laplacian MEDI
        set(h.phaseUnwrap.checkbox.excludeMask, 'Enable', 'off', 'Value', 0);
        set(h.phaseUnwrap.edit.excludeMask,     'Enable', 'off');
        
    case methodUnwrapName{2}    % Laplacian STI suite
        set(h.phaseUnwrap.checkbox.excludeMask, 'Enable', 'off', 'Value', 0);
        set(h.phaseUnwrap.edit.excludeMask,     'Enable', 'off');

    case methodUnwrapName{3}    % 3D best path
        set(h.phaseUnwrap.checkbox.excludeMask, 'Enable', 'on');

    case methodUnwrapName{4}    % region growing medi
        set(h.phaseUnwrap.checkbox.excludeMask, 'Enable', 'on');

    case methodUnwrapName{5}    % graphcut
        set(h.phaseUnwrap.checkbox.excludeMask, 'Enable', 'on');
        
    case methodUnwrapName{6}    % segue
        set(h.phaseUnwrap.checkbox.excludeMask, 'Enable', 'on');
    
    % in the future, add new method here
end

end

% echo phase combination method specific panel setting
function popupPhaseCombine_Callback(source,eventdata,h)

sepia_universal_variables;

% get selected background removal method
method = source.String{source.Value,1} ;

% switch on target panel
switch method
    case methodEchoCombineName{1}   % optimal weighted combination
        set(h.phaseUnwrap.checkbox.saveEchoPhase,   'Enable', 'on', 'Value', 0);
        set(h.phaseUnwrap.checkbox.eddyCorrect,     'Enable', 'on');
        
    case methodEchoCombineName{2}   % MEDI nonlinear fit
        set(h.phaseUnwrap.checkbox.saveEchoPhase,   'Enable', 'off', 'Value', 0);
        set(h.phaseUnwrap.checkbox.eddyCorrect,     'Enable', 'on');
        
    case methodEchoCombineName{3}   % MEDI nonlinear fit with bipolar readout
        set(h.phaseUnwrap.checkbox.saveEchoPhase,   'Enable', 'off', 'Value', 0);
        set(h.phaseUnwrap.checkbox.eddyCorrect,     'Enable', 'off', 'Value', 0);
    
end

end