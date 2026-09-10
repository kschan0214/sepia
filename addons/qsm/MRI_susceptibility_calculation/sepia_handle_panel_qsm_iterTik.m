%% h = sepia_handle_panel_qsm_iterTik(hParent,h,position)
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
% Description: This GUI function creates a panel for iterTik method
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 4 August 2020
% Date modified: 
%
%
function h = sepia_handle_panel_qsm_iterTik(hParent,h,position)

%% set default values
defaultTol          = 0.03;
defaultThreshold	= 2/3;
defaultLambda       = 0.05;

menuSolver       = {'Truncated kspace division','Direct Tikhonov','Iterative Tikhonov'};

%% Tooltips
tooltip.qsm.iterTik.tol            = 'Convergence limit, change rate in the solution';
tooltip.qsm.iterTik.threshold      = 'Conjugate gradient stopping tolerance';
tooltip.qsm.iterTik.lambda         = 'Regularisation value';

%% layout of the panel
nrow        = 4;
rspacing    = 0.03;
ncol        = 2;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent handle of iterTik panel children

h.qsm.panel.iterTik = uipanel(hParent,...
        'Title','Truncated K-space division + direct/iterative Tikhonov regularisation',...
        'position',position,...
        'backgroundcolor',get(h.fig,'color'),'Visible','off');

%% Children of iterTik panel
    
    panelParent = h.qsm.panel.iterTik;

    % width of each element in a functional column, in normalised unit
    wratio = 0.5;
    
    % col 1, row 1
    % text|popup field pair: solver
    [h.qsm.iterTik.text.solver,h.qsm.iterTik.popup.solver] = sepia_construct_text_popup(...
        panelParent,'Solver:',                  menuSolver,         [left(1) bottom(1) width height], wratio);
    
    % col 1, row 2
    % text|edit field pair: maximum iterations
    [h.qsm.iterTik.text.threshold,h.qsm.iterTik.edit.threshold] = sepia_construct_text_edit(...
        panelParent,'K-space threshold:',       defaultThreshold,  	[left(1) bottom(2) width height], wratio);
    
    % col 1, row 3
    % text|edit field pair: maximum iterations
    [h.qsm.iterTik.text.lambda,h.qsm.iterTik.edit.lambda] = sepia_construct_text_edit(...
        panelParent,'Lambda:',                  defaultLambda,      [left(1) bottom(3) width height], wratio);
    set(h.qsm.iterTik.edit.lambda, 'Enable', 'off');
    
    % col 1, row 4
    % text|edit field pair: tolerance
    [h.qsm.iterTik.text.tol,h.qsm.iterTik.edit.tol] = sepia_construct_text_edit(...
        panelParent,'CG Tolerance:',          	defaultTol,         [left(1) bottom(4) width height], wratio);
    set(h.qsm.iterTik.edit.tol, 'Enable', 'off');
    
    
%% set tooltips
set(h.qsm.iterTik.text.tol,               'Tooltip',tooltip.qsm.iterTik.tol);
set(h.qsm.iterTik.text.threshold,         'Tooltip',tooltip.qsm.iterTik.threshold);
set(h.qsm.iterTik.text.lambda,            'Tooltip',tooltip.qsm.iterTik.lambda);

%% set callbacks
set(h.qsm.iterTik.edit.lambda,          'Callback', {@EditInputMinMax_Callback,defaultLambda,   0,0});
set(h.qsm.iterTik.edit.threshold,       'Callback', {@EditInputMinMax_Callback,defaultThreshold,0,0});
set(h.qsm.iterTik.edit.tol,             'Callback', {@EditInputMinMax_Callback,defaultTol,      0,0});
set(h.qsm.iterTik.popup.solver,         'Callback', {@switch_solver_callback,h,menuSolver});

end

% Callback for switching solver
function switch_solver_callback(source,eventdata,h,menuSolver)

solver = source.String{source.Value};

switch solver
    case menuSolver{1}
        set(h.qsm.iterTik.edit.threshold,	'Enable', 'on');
        set(h.qsm.iterTik.edit.lambda,      'Enable', 'off');
        set(h.qsm.iterTik.edit.tol,         'Enable', 'off');

    case menuSolver{2}
        set(h.qsm.iterTik.edit.threshold,  	'Enable', 'off');
        set(h.qsm.iterTik.edit.lambda,      'Enable', 'on');
        set(h.qsm.iterTik.edit.tol,         'Enable', 'off');

    case menuSolver{3}
        set(h.qsm.iterTik.edit.threshold,	'Enable', 'off');
        set(h.qsm.iterTik.edit.lambda,     	'Enable', 'on');
        set(h.qsm.iterTik.edit.tol,      	'Enable', 'on');
end

% Truncated kspace division / Direct Tikhonov apply the mask only after a
% closed-form/direct k-space inversion, so two-pass masking has no effect
% on their output by construction - warn if it's currently enabled (this
% panel is only visible when the top-level QSM method is 'MRI Suscep.
% Calc.' itself, so no need to re-check that here)
if isfield(h.qsm,'popup') && isfield(h.qsm.popup,'twopass')
    twopass = h.qsm.popup.twopass.String{h.qsm.popup.twopass.Value};
    if ~strcmpi(twopass,'None') && ~qsm_method_uses_mask_in_inversion('MRI Suscep. Calc.', solver)
        warndlg(sprintf(['The selected solver (MRI Suscep. Calc. / %s) applies the mask only ', ...
            'after a closed-form/direct k-space inversion, so two-pass masking (%s) will have ', ...
            'no effect on the result.'], solver, twopass), ...
            'Two-pass masking has no effect with this method');
    end
end

end