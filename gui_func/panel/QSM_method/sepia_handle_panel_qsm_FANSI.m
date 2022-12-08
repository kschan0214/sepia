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
% Description: This GUI function creates a panel for FANSI method
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 1 June 2018
% Date modified: 4 April 2020 (v0.8.0)
% Date modified: 20 Feb 2022 (v1.0)
%
%
function h = sepia_handle_panel_qsm_FANSI(hParent,h,position)

%% set default values
% defaultTol      = 1;  
% defaultLambda   = 3e-5;
% defaultMu       = 5e-5; 
% defaultMaxIter  = 50;  
% defaultBeta = 1;
% defaultMuh  = 1;

defaultMu2      = 1;
defaultLambda   = 2e-4;              % update v3
defaultTol      = 0.1;               % update v3
defaultMu       = defaultLambda*100; % update v3
defaultMaxIter  = 150;               % update v3
defaultBeta = 150;                   % update v3
defaultMuh  = 3;                     % update v3

menuSolver       = {'Non-linear','Linear'};
menuConstraints  = {'TGV','TV'};
menuGradientMode = {'Vector field','L1 norm','L2 norm','None'};

%% Tooltips
tooltip.qsm.fansi.tol            = 'Convergence limit, change rate in the solution, adpated for FANSI v3';
tooltip.qsm.fansi.MaxIter        = 'Maximum iterations allowed, adpated for FANSI v3';
tooltip.qsm.fansi.lambda         = 'i.e. alpha1, adpated for FANSI v3';
tooltip.qsm.fansi.mu             = 'i.e. mu1, adpated for FANSI v3; recommended value=100*alpha0';
tooltip.qsm.fansi.mu2            = 'i.e. mu2; recommended value=1';
tooltip.qsm.fansi.isWeakHarmonic = 'Check to us weak-harmonic regularization algorithms';
tooltip.qsm.fansi.beta           = 'i.e. beta, adpated for FANSI v3';
tooltip.qsm.fansi.muh            = 'i.e. muh, adpated for FANSI v3; recommended value=beta/50';

%% layout of the panel
nrow        = 4;
rspacing    = 0.03;
ncol        = 3;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent handle of FANSI panel children

h.qsm.panel.FANSI = uipanel(hParent,...
        'Title','FAst Nonlinear Susceptibility Inversion (FANSI)',...
        'position',position,...
        'backgroundcolor',get(h.fig,'color'),'Visible','off');

%% Children of FANSI panel
    
    panelParent = h.qsm.panel.FANSI;

    % width of each element in a functional column, in normalised unit
    wratio = 0.65;
    
    % col 1, row 1
    % text|edit field pair: tolerance
    [h.qsm.FANSI.text.tol,h.qsm.FANSI.edit.tol] = sepia_construct_text_edit(...
        panelParent,'Tolerance:',               defaultTol,         [left(1) bottom(1) width height], wratio);
    
    % col 1, row 2
    % text|edit field pair: maximum iterations
    [h.qsm.FANSI.text.maxIter,h.qsm.FANSI.edit.maxIter] = sepia_construct_text_edit(...
        panelParent,'Max. iterations:',         defaultMaxIter,     [left(1) bottom(2) width height], wratio);
    
    % col 1, row 3
    % text|popup field pair: solver
    [h.qsm.FANSI.text.solver,h.qsm.FANSI.popup.solver] = sepia_construct_text_popup(...
        panelParent,'Solver:',                  menuSolver,         [left(1) bottom(3) width height], wratio);
    
    % col 1, row 4
    % text|popup field pair: constraint
    [h.qsm.FANSI.text.constraints,h.qsm.FANSI.popup.constraints] = sepia_construct_text_popup(...
        panelParent,'Constraint:',              menuConstraints,	[left(1) bottom(4) width height], wratio);
    
    % col 2, row 1
    % text|edit field pair: fidelity consistency
    [h.qsm.FANSI.text.mu2,h.qsm.FANSI.edit.mu2] = sepia_construct_text_edit(...
        panelParent,'Fidelity consistency:',	defaultMu2,         [left(2) bottom(1) width height], wratio);
    
    % col 2, row 2
    % text|edit field pair: gradient L1 penalty
    [h.qsm.FANSI.text.lambda,h.qsm.FANSI.edit.lambda] = sepia_construct_text_edit(...
        panelParent,'Gradient L1 penalty:',     defaultLambda,      [left(2) bottom(2) width height], wratio);
    
    % col 2, row 3
    % text|edit field pair: gradient consistency
    [h.qsm.FANSI.text.mu,h.qsm.FANSI.edit.mu] = sepia_construct_text_edit(...
        panelParent,'Gradient consistency:',    defaultMu,          [left(2) bottom(3) width height], wratio);

    % col 2, row 4
    % text|popup field pair: gradient mode
    [h.qsm.FANSI.text.gradientMode,h.qsm.FANSI.popup.gradientMode] = sepia_construct_text_popup(...
        panelParent,'Gradient mode:',           menuGradientMode,	[left(2) bottom(4) width height], wratio);
    

    % col 3, row 1
    % weak harmonic field 
    h.qsm.FANSI.checkbox.isWeakHarmonic = uicontrol('Parent',h.qsm.panel.FANSI ,...
        'Style','checkbox','String','Weak-Harmonic Regularisation',...
        'units','normalized','Position',[left(3) bottom(1) width height],...
        'backgroundcolor',get(h.fig,'color'));
    
    % col 3, row 2
    % text|edit field pair: harmonic constrain weight
    [h.qsm.FANSI.text.beta,h.qsm.FANSI.edit.beta] = sepia_construct_text_edit(...
        panelParent,'Harmonic constraint:',     defaultBeta,        [left(3) bottom(2) width height], wratio);
    set(h.qsm.FANSI.edit.beta, 'Enable', 'off');
    
    % col 3, row 3
    % text|edit field pair: harmonic consistency weight
    [h.qsm.FANSI.text.muh,h.qsm.FANSI.edit.muh] = sepia_construct_text_edit(...
        panelParent,'Harmonic consistency:',    defaultMuh,         [left(3) bottom(3) width height], wratio);
    set(h.qsm.FANSI.edit.muh, 'Enable', 'off');

    % col 3, row 4
    % checkbox field pair: isGPU
    h.qsm.FANSI.checkbox.isGPU = uicontrol('Parent',h.qsm.panel.FANSI,'Style','checkbox',...
        'String','Enable GPU',...
        'units','normalized','position',[left(3) bottom(4) width height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'));
    
    
%% set tooltips
set(h.qsm.FANSI.text.tol,               'Tooltip',tooltip.qsm.fansi.tol);
set(h.qsm.FANSI.text.maxIter,           'Tooltip',tooltip.qsm.fansi.MaxIter);
set(h.qsm.FANSI.text.lambda,            'Tooltip',tooltip.qsm.fansi.lambda);
set(h.qsm.FANSI.checkbox.isWeakHarmonic,'Tooltip',tooltip.qsm.fansi.isWeakHarmonic);
set(h.qsm.FANSI.text.mu,                'Tooltip',tooltip.qsm.fansi.mu);
set(h.qsm.FANSI.text.beta,              'Tooltip',tooltip.qsm.fansi.beta);
set(h.qsm.FANSI.text.muh,               'Tooltip',tooltip.qsm.fansi.muh);
set(h.qsm.FANSI.text.mu2,               'Tooltip',tooltip.qsm.fansi.mu2);

%% set callbacks
% set(h.qsm.FANSI.edit.lambda,	'Callback', {@EditInputMinMax_Callback,defaultLambda,   0,0});
set(h.qsm.FANSI.edit.lambda,	'Callback', {@EditInputMinMax_lambda2mu_Callback,defaultLambda,   0,0,[],h});
set(h.qsm.FANSI.edit.maxIter,	'Callback', {@EditInputMinMax_Callback,defaultMaxIter,  1,1});
set(h.qsm.FANSI.edit.mu,     	'Callback', {@EditInputMinMax_Callback,defaultMu,       0,0});
set(h.qsm.FANSI.edit.mu2,     	'Callback', {@EditInputMinMax_Callback,defaultMu2,      0,0});
set(h.qsm.FANSI.edit.tol,    	'Callback', {@EditInputMinMax_Callback,defaultTol,      0,0});
set(h.qsm.FANSI.edit.beta,     	'Callback', {@EditInputMinMax_Callback,defaultBeta,     0,0});
set(h.qsm.FANSI.edit.muh,     	'Callback', {@EditInputMinMax_Callback,defaultMuh,      0,0});

set(h.qsm.FANSI.checkbox.isWeakHarmonic,	'Callback', {@CheckboxEditPair_Callback,{h.qsm.FANSI.edit.beta,h.qsm.FANSI.edit.muh},1});

end

function EditInputMinMax_lambda2mu_Callback(source,eventdata,defaultValue,isIntegerInput,lb,ub,h)
% set the min/max input allowed in edit fields

    % check input is numeric
    if isnan(str2double(source.String)) || ~isreal(str2double(source.String))
        warndlg('Please enter a valid real number');
        source.String = num2str(defaultValue);
    end

    % check minimum
    if str2double(source.String)<lb
        source.String = num2str(lb);
        warndlg(['The minimium value allowed is ' num2str(lb)]);
    end
    
    % if maximum is set then check maximum
    if ~isempty(ub)
        if str2double(source.String)>ub
            source.String = num2str(ub);
            warndlg(['The maximum value allowed is ' num2str(ub)]);
        end
    end
    
    % make sure the input is interger for some fields
    if isIntegerInput
        source.String = num2str(round(str2double(source.String)));
    end
    
    set_non_nan_value(h.qsm.FANSI.edit.mu,'String',num2str(100*str2double(source.String)))
    
end