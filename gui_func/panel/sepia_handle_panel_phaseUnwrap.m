%% h = sepia_handle_panel_phaseUnwrap(hParent,h,position)
%
% Input
% --------------
% hParent       : parent handle of this panel
% hFig          : handle of the GUI
% h             : global structure contains all handles
% position      : position of this panel
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
%
%
function h = sepia_handle_panel_phaseUnwrap(hParent,h,position)

% set up method name displayed on GUI
methodUnwrapName        = {'Laplacian','Laplacian STI suite','3D best path','Region growing','Graphcut','SEGUE'};
methodEchoCombineName   = {'Optimum weights','MEDI nonlinear fit','MEDI nonlinear fit (Bipolar)'};
methodExcludedName      = {'Weighting map','Brain mask'};

% set default value
defaultThreshold = 0.5;

% define maximum level of options and spacing between options
nlevel = 5;
spacing = 0.02;
height = (1-(nlevel+1)*spacing)/nlevel;
button = (height+spacing:height+spacing:(height+spacing)*nlevel) - height;

% Parent handle of phase unwrapping panel
h.StepsPanel.phaseUnwrap = uipanel(hParent,'Title','Total field recovery and phase unwrapping',...
    'position',[position(1) position(2) 0.95 0.2]);
    
    % Temporo-spatial unwrapping methods
    h.phaseUnwrap.text.phaseCombMethod = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,...
        'Style','text','String','Echo phase combination:',...
        'units','normalized','position',[0.01 button(nlevel) 0.3 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Select method to compute total field map');
    h.phaseUnwrap.popup.phaseCombMethod = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,...
        'Style','popup',...
        'String',methodEchoCombineName,...
        'units','normalized','position',[0.31 button(nlevel) 0.4 height]); 
    
    % phase unwrapping method
    h.phaseUnwrap.text.phaseUnwrap = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,...
        'Style','text','String','Phase unwrapping:',...
        'units','normalized','position',[0.01 button(nlevel-1) 0.3 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Select phase unwrapping algorithm');
    h.phaseUnwrap.popup.phaseUnwrap = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,...
        'Style','popup',...
        'String',methodUnwrapName,...
        'units','normalized','position',[0.31 button(nlevel-1) 0.4 height]); 

    % eddy current correction for bipolar readout
    h.phaseUnwrap.checkbox.eddyCorrect = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,...
        'Style','checkbox','String','Bipolar readout eddy current correction',...
        'units','normalized','Position',[0.01 button(nlevel-2) 1 height],...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Correct inconsistent phase between odd and even number echoes');
    
    % exclusion of unreliable voxels
    h.phaseUnwrap.checkbox.excludeMask = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,...
        'Style','checkbox','String','Exclude unreliable voxels, threshold:',...
        'units','normalized','position',[0.01 button(nlevel-3) 0.3 height],...
        'backgroundcolor',get(h.fig,'color'),...
        'Enable','off',...
        'tooltip','Apply threshold to exclude unreliable voxels based on relative residual of mono-exponential decay model with a single frequency shift (only available for region growing based methods)');
    h.phaseUnwrap.edit.excludeMask = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,...
        'Style','edit',...
        'String',num2str(defaultThreshold),...
        'units','normalized','position',[0.31 button(nlevel-3) 0.04 height],...
        'backgroundcolor','white',...
        'Enable','off',...
        'tooltip','Higher value means accepting larger error between the data fitting and measurement');
    % excluding method
    h.phaseUnwrap.text.excludeMethod = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,...
        'Style','text','String','in ',...
        'units','normalized','position',[0.36 button(nlevel-3) 0.03 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'));
    h.phaseUnwrap.popup.excludeMethod = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,...
        'Style','popup',...
        'String',methodExcludedName,...
        'Enable','off',...
        'units','normalized','position',[0.40 button(nlevel-3) 0.3 height]); 
    
    % save unwrapped echo phase option
    h.phaseUnwrap.checkbox.saveEchoPhase = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,...
        'Style','checkbox','String','Save unwrapped echo phase',...
        'units','normalized','position',[0.01 button(nlevel-4) 0.3 height],...
        'backgroundcolor',get(h.fig,'color'),...
        'Enable','on',...
        'tooltip','If total field is computed using optimal weighted combination with multi echo data, it is possible to save the unwrapped phase for each echo.');

%% set callback functions
set(h.phaseUnwrap.checkbox.excludeMask,	'Callback', {@CheckboxEditPair_Callback,{h.phaseUnwrap.edit.excludeMask,h.phaseUnwrap.popup.excludeMethod},1});
set(h.phaseUnwrap.edit.excludeMask, 	'Callback', {@EditInputMinMax_Callback,defaultThreshold,0,0,1});
set(h.phaseUnwrap.popup.phaseUnwrap,    'Callback', {@popupPhaseUnwrap_Callback,h});
set(h.phaseUnwrap.popup.phaseCombMethod,'Callback', {@popupPhaseCombine_Callback,h});
        
end

%% Callback functions
% Phase unwrapping method
function popupPhaseUnwrap_Callback(source,eventdata,h)

% get selected background removal method
method = source.String{source.Value,1} ;

% switch on target panel
switch method
    case 'Laplacian'
        set(h.phaseUnwrap.checkbox.excludeMask, 'Enable', 'off', 'Value', 0);
        set(h.phaseUnwrap.edit.excludeMask,     'Enable', 'off');
        
    case 'Laplacian STI suite'
        set(h.phaseUnwrap.checkbox.excludeMask, 'Enable', 'off', 'Value', 0);
        set(h.phaseUnwrap.edit.excludeMask,     'Enable', 'off');

    case '3D best path'
        set(h.phaseUnwrap.checkbox.excludeMask, 'Enable', 'on');

    case 'Region growing'
        set(h.phaseUnwrap.checkbox.excludeMask, 'Enable', 'on');

    case 'Graphcut'
        set(h.phaseUnwrap.checkbox.excludeMask, 'Enable', 'on');
        
    case 'SEGUE'
        set(h.phaseUnwrap.checkbox.excludeMask, 'Enable', 'on');
    
    % in the future, add new method here
end

end

% echo phase combination method
function popupPhaseCombine_Callback(source,eventdata,h)

% get selected background removal method
method = source.String{source.Value,1} ;

% switch on target panel
switch method
    case 'Optimum weights'
        set(h.phaseUnwrap.checkbox.saveEchoPhase,   'Enable', 'on', 'Value', 0);
        set(h.phaseUnwrap.checkbox.eddyCorrect,     'Enable', 'on');
        
    case 'MEDI nonlinear fit'
        set(h.phaseUnwrap.checkbox.saveEchoPhase,   'Enable', 'off', 'Value', 0);
        set(h.phaseUnwrap.checkbox.eddyCorrect,     'Enable', 'on');
        
    case 'MEDI nonlinear fit (Bipolar)'
        set(h.phaseUnwrap.checkbox.saveEchoPhase,   'Enable', 'off', 'Value', 0);
        set(h.phaseUnwrap.checkbox.eddyCorrect,     'Enable', 'off', 'Value', 0);
    
end

end