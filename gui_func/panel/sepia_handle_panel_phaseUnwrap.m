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
% Date modified: 28 June 2020 (v0.8.0)
%
%
function h = sepia_handle_panel_phaseUnwrap(hParent,h,position)

% load universal variable for methods' names
sepia_universal_variables;

% set default value
defaultThreshold = 0.5;

% tooltips
tooltip.unwrap.panel.method         = 'Select a method to combine field maps from multi-echo data';
tooltip.unwrap.panel.unwrap         = 'Select a phase unwrapping algorithm for spatial unwrapping';
tooltip.unwrap.panel.eddy           = 'Correct phase offset between odd and even number echoes if bipolar readout was used in acquisition';
tooltip.unwrap.panel.exclude        = ['Apply threshold on relative residual to exclude unreliable voxels based on ',...
                                        'mono-exponential decay model (only available for non-Laplacian methods)'];
tooltip.unwrap.panel.exclude_edit	= ['Higher value means accepting larger error between the data fitting and measurement'];
tooltip.unwrap.panel.saveUnwrap     = ['If total field is computed using optimal weights combination with multi-echo data, ',...
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
    wratio = 0.5;
    
    % row 1
    % Temporo-spatial unwrapping methods, 'text|popup' functional column,
    [h.phaseUnwrap.text.phaseCombMethod,h.phaseUnwrap.popup.phaseCombMethod] = sepia_construct_text_popup(...
        h.StepsPanel.phaseUnwrap,'Echo phase combination:', methodEchoCombineName, [left(1) bottom(1) width height], wratio);

    % row 2
    % phase unwrapping method, 'text|popup' 
    [h.phaseUnwrap.text.phaseUnwrap,h.phaseUnwrap.popup.phaseUnwrap] = sepia_construct_text_popup(...
        h.StepsPanel.phaseUnwrap,'Phase unwrapping:', methodUnwrapName, [left(1) bottom(2) width height], wratio);
    
    % row 3
    % eddy current correction for bipolar readout, 'checkbox' functional
    h.phaseUnwrap.checkbox.eddyCorrect = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,...
        'Style','checkbox','String','Bipolar readout correction',...
        'units','normalized','Position',[0.01 bottom(3) 1 height],...
        'backgroundcolor',get(h.fig,'color'));
    
    % exclusion of unreliable voxels, 'checkbox|field|text|popup', 4th row
    h.phaseUnwrap.checkbox.excludeMask = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,...
        'Style','checkbox','String','Exclude voxels using residual, threshold:',...
        'units','normalized','position',[0.01 bottom(4) 0.3 height],...
        'backgroundcolor',get(h.fig,'color'),...
        'Enable','off');
    h.phaseUnwrap.edit.excludeMask = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,...
        'Style','edit',...
        'String',num2str(defaultThreshold),...
        'units','normalized','position',[0.31 bottom(4) 0.04 height],...
        'backgroundcolor','white',...
        'Enable','off');
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
        'Enable','on');

%% set tooltips
set(h.phaseUnwrap.text.phaseCombMethod,     'Tooltip',tooltip.unwrap.panel.method);
set(h.phaseUnwrap.text.phaseUnwrap,         'Tooltip',tooltip.unwrap.panel.unwrap);
set(h.phaseUnwrap.checkbox.eddyCorrect,     'Tooltip',tooltip.unwrap.panel.eddy);
set(h.phaseUnwrap.checkbox.excludeMask,     'Tooltip',tooltip.unwrap.panel.exclude);
set(h.phaseUnwrap.edit.excludeMask,         'Tooltip',tooltip.unwrap.panel.exclude_edit);
set(h.phaseUnwrap.checkbox.saveEchoPhase,	'Tooltip',tooltip.unwrap.panel.saveUnwrap);

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

% Reset the option 
set(h.phaseUnwrap.checkbox.excludeMask, 'Enable', 'off', 'Value', 0);
set(h.phaseUnwrap.edit.excludeMask,     'Enable', 'off');
set(h.phaseUnwrap.popup.excludeMethod,  'Enable', 'off');
% method the user chosen will affect if exclusion method can be used or not 
for k = 1:length(methodUnwrapName)
    if strcmpi(method,methodUnwrapName{k})
        set(h.phaseUnwrap.checkbox.excludeMask, 'Enable', gui_unwrap_exclusion{k});
    end
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