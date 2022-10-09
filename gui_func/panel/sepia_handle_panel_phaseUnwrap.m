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
% Date modified: 12 June 2021 (v1.0)
%
%
function h = sepia_handle_panel_phaseUnwrap(hParent,h,position)

% load universal variable for methods' names
sepia_universal_variables;

% set default value
defaultThreshold = 0.5;

% tooltips
tooltip.unwrap.panel.method         = 'Select a method to combine field maps from multi-echo data';

%% layout of the panel
nrow        = 4;
rspacing    = 0.01;
ncol        = 2;
cspacing    = 0.01;
[~,~,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

% Set parent of phase unwrapping panel
h.StepsPanel.phaseUnwrap = uipanel(hParent,...
    'Title','Total field recovery and phase unwrapping',...
    'fontweight', 'bold',...
    'position',[position(1) position(2) 0.95 0.2]);
%% design of this panel

    
    % width of each element in a functional column, in normalised unit of
    % the functional column width
    height = 0.1;   % row height
    wratio = 0.5;   % colume width
    
    % row 1
    % Temporo-spatial unwrapping methods, 'text|popup' functional column,
    % col 1
    % text|popup pair: select method
    [h.phaseUnwrap.text.phaseCombMethod,h.phaseUnwrap.popup.phaseCombMethod] = sepia_construct_text_popup(...
        h.StepsPanel.phaseUnwrap,'Echo phase combination:', methodEchoCombineName, [left(1) 0.85 width height], wratio);

%% create control panel

% define position and size of all method panels
position_child = [0.01 0.04 0.95 0.75];

% construct all method panels
for k = 1:length(function_EchoCombine_method_panel)
    h = feval(function_EchoCombine_method_panel{k},h.StepsPanel.phaseUnwrap,h,position_child);
end
 
%% set tooltips
set(h.phaseUnwrap.text.phaseCombMethod, 'Tooltip',tooltip.unwrap.panel.method);
% 
%% set callback functions
set(h.phaseUnwrap.popup.phaseCombMethod, 'Callback', {@PopupEchoCombine_Callback,h});
        
end

%% Callback functions
% display corresponding QSM method's panel
function PopupEchoCombine_Callback(source,eventdata,h)

% needs 'methodQSMName' here
sepia_universal_variables;

% get selected QSM method
method = source.String{source.Value,1} ;

% switch off all panels
fields = fieldnames(h.phaseUnwrap.panel);
for kf = 1:length(fields)
    set(h.phaseUnwrap.panel.(fields{kf}),   'Visible','off');
end

% switch on only target panel
for k = 1:length(methodEchoCombineName)
    if strcmpi(method,methodEchoCombineName{k})
        set(h.phaseUnwrap.panel.(fields{k}), 'Visible','on');
        break
    end
end

end
