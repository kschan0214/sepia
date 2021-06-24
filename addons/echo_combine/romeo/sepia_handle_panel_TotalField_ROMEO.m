%% h = sepia_handle_panel_TotalField_ROMEO(hParent,h,position)
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
% Description: This GUI function creates a panel for ROMEO method
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 12 June 2021 (v1.0)
% Date modified: 
%
%
function h = sepia_handle_panel_TotalField_ROMEO(hParent,h,position)
% set up method name displayed on GUI
sepia_universal_variables;

%% set default values
defaultThreshold = 0.5;
%defaultMcpc3dsSelection = 1;

%% Tooltips
% tooltips
tooltip.unwrap.panel.unwrap         = 'Select a phase unwrapping algorithm for spatial unwrapping';
tooltip.unwrap.panel.exclude        = ['Apply threshold on relative residual to exclude unreliable voxels based on ',...
                                        'mono-exponential decay model (only available for non-Laplacian methods)'];
tooltip.unwrap.panel.exclude_edit	= ['Higher value means accepting larger error between the data fitting and measurement'];

%% layout of the panel
nrow        = 3;
rspacing    = 0.01;
ncol        = 2;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent handle of children panel 

h.phaseUnwrap.panel.ROMEOTotalField = uipanel(hParent,...
    'Title','ROMEO total field calculation',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','off');

%% Children of panel
    
    panelParent = h.phaseUnwrap.panel.ROMEOTotalField;

    % width of each element in a functional column, in normalised unit
    wratio = 0.5;
    
    krow = 1;
    % row 1, left
    % MCPC-3D-S phase offset correction popup
    methodOffsetCorrect = {'Off', 'On', 'Bipolar (>= 3 echoes)'};
    [h.phaseUnwrap.ROMEOTotalField.text.offsetCorrect,h.phaseUnwrap.ROMEOTotalField.popup.offsetCorrect] = sepia_construct_text_popup(...
        panelParent,'MCPC-3D-S phase offset correction:', methodOffsetCorrect, [left(1) bottom(krow) width height], wratio);
    
    % row 1, right
    % eddy current correction for bipolar readout, 'checkbox' functional
    h.phaseUnwrap.ROMEOTotalField.checkbox.eddyCorrect = uicontrol('Parent',panelParent ,...
        'Style','checkbox','String','SEPIA Bipolar readout correction (alternative to MPCP-3D-S bipolar)',...
        'units','normalized','Position',[left(2) bottom(krow) width height],...
        'backgroundcolor',get(h.fig,'color'));
    
    krow = krow + 1;
    % row 2, left
    % Mask used for ROMEO
    methodMask = {'SEPIA mask', 'ROMEO mask', 'No Mask'};
    [h.phaseUnwrap.ROMEOTotalField.text.mask,h.phaseUnwrap.ROMEOTotalField.popup.mask] = sepia_construct_text_popup(...
        panelParent,'Mask for unwrapping:', methodMask, [left(1) bottom(krow) width height], wratio);
    
    % row 2, right
    % save unwrapped echo phase option, 'checkbox'
    h.phaseUnwrap.ROMEOTotalField.checkbox.saveEchoPhase = uicontrol('Parent',panelParent ,...
        'Style','checkbox','String','Save unwrapped echo phase',...
        'units','normalized','position',[left(2) bottom(krow) width height],...
        'backgroundcolor',get(h.fig,'color'),...
        'Enable','on');
    
    % row 3
    krow = krow + 1;
    % exclusion of unreliable voxels, 'checkbox|field|text|popup'
    h.phaseUnwrap.ROMEOTotalField.checkbox.excludeMask = uicontrol('Parent',panelParent ,...
        'Style','checkbox','String','Exclude voxels using residual, threshold:',...
        'units','normalized','position',[left(1) bottom(krow) 0.4 height],...
        'backgroundcolor',get(h.fig,'color'),...
        'Enable','off');
    h.phaseUnwrap.ROMEOTotalField.edit.excludeMask = uicontrol('Parent',panelParent ,...
        'Style','edit',...
        'String',num2str(defaultThreshold),...
        'units','normalized','position',[left(1)+0.4 bottom(krow) 0.04 height],...
        'backgroundcolor','white',...
        'Enable','off');
    % excluding method
    h.phaseUnwrap.ROMEOTotalField.text.excludeMethod = uicontrol('Parent',panelParent ,...
        'Style','text','String','and apply in ',...
        'units','normalized','position',[left(1)+0.46 bottom(krow) 0.1 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'));
    h.phaseUnwrap.ROMEOTotalField.popup.excludeMethod = uicontrol('Parent',panelParent ,...
        'Style','popup',...
        'String',methodExcludedName,...
        'Enable','off',...
        'units','normalized','position',[left(1)+0.56 bottom(krow) 0.23 height]); 
    

%% set tooltips
%set(h.phaseUnwrap.ROMEOTotalField.text.phaseUnwrap,          'Tooltip',tooltip.unwrap.panel.unwrap);
set(h.phaseUnwrap.ROMEOTotalField.checkbox.excludeMask,      'Tooltip',tooltip.unwrap.panel.exclude);
set(h.phaseUnwrap.ROMEOTotalField.edit.excludeMask,          'Tooltip',tooltip.unwrap.panel.exclude_edit);

%% set callbacks
set(h.phaseUnwrap.ROMEOTotalField.checkbox.excludeMask,  'Callback', {@CheckboxEditPair_Callback,{h.phaseUnwrap.ROMEOTotalField.edit.excludeMask,h.phaseUnwrap.ROMEOTotalField.popup.excludeMethod},1});
set(h.phaseUnwrap.ROMEOTotalField.edit.excludeMask,      'Callback', {@EditInputMinMax_Callback,defaultThreshold,0,0,1});
%set(h.phaseUnwrap.ROMEOTotalField.popup.phaseUnwrap,     'Callback', {@popupPhaseUnwrap_Callback,h});

% Enable exlude options 
set(h.phaseUnwrap.ROMEOTotalField.checkbox.excludeMask, 'Enable', 'on', 'Value', 0);
set(h.phaseUnwrap.ROMEOTotalField.edit.excludeMask,     'Enable', 'on');
set(h.phaseUnwrap.ROMEOTotalField.popup.excludeMethod,  'Enable', 'on');

end

%% Callback functions
% Phase unwrapping method specific panel setting
function popupPhaseUnwrap_Callback(source,eventdata,h)

sepia_universal_variables;

% get selected background removal method
method = source.String{source.Value,1} ;

% Reset the option 
set(h.phaseUnwrap.optimumWeights.checkbox.excludeMask, 'Enable', 'off', 'Value', 0);
set(h.phaseUnwrap.optimumWeights.edit.excludeMask,     'Enable', 'off');
set(h.phaseUnwrap.optimumWeights.popup.excludeMethod,  'Enable', 'off');
% method the user chosen will affect if exclusion method can be used or not 
for k = 1:length(methodUnwrapName)
    if strcmpi(method,methodUnwrapName{k})
        set(h.phaseUnwrap.optimumWeights.checkbox.excludeMask, 'Enable', gui_unwrap_exclusion{k});
    end
end

end
