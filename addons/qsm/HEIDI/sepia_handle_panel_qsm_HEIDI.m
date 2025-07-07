%% h = sepia_handle_panel_qsm_HEIDI(hParent,h,position)
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
% Description: This GUI function creates a panel for HEIDI method
%
% Kwok-shing Chan @ MGH
% kchan2@mgh.harvard.edu
% Date created: 14 June 2025
% Date modified: 
%
%
function h = sepia_handle_panel_qsm_HEIDI(hParent,h,position)

%% set default values   
%%%%%%%% This are the default value on the GUI
defaultTol                      = 1e-5;              
defaultMaxIter                  = 400;               
defaultOffsetUseBool            = true;         
defaultIsFourierDomainFormula   = false;   
defaultProcProcConeThreshold    = 0.1;
defaultProcProcConeTol          = eps;
defaultProcProcConeTolEnergy    = eps;
defaultResidualWeighting        = 0.2 / 9.4; % will be scale by field strength

% full lists, needs further testing, current not available
% menuTikhonov        = {'Default','Partial Gradient weighting','Laplacian'};
% menuSolver          = {'Default','InverseFiltering','SpatialDomainTV'};
% menuDipoleFilter    = {'Default','truncSingularValues'};
menuTikhonov        = {'Default'};
menuSolver          = {'Default'};
menuDipoleFilter    = {'Default'};

%% Tooltips
%%%%%%%% (Optional) Tooltips to elaborate algorithm parameter usage
tooltip.qsm.heidi.tol                       = 'Convergence limit';
tooltip.qsm.heidi.MaxIter                   = 'Maximum iterations allowed';
tooltip.qsm.heidi.isOffsetUse               = 'Automatic offset determination';
tooltip.qsm.heidi.isFourierDomainFormula    = '';
tooltip.qsm.heidi.Tikhonov                  = 'Performing dipole inversion with Tikhonov regularisation';
tooltip.qsm.heidi.solver                    = '';
tooltip.qsm.heidi.DipoleFilter              = '';
tooltip.qsm.heidi.ResidualWeighting         = 'This value will be scaled by field strength, i.e. ResidualWeighting*B0';
tooltip.qsm.heidi.ProcProcConeThreshold     = 'Threshold for post-processing of cone';
tooltip.qsm.heidi.ProcProcConeTol           = 'Tolerance for post-processing of cone';
tooltip.qsm.heidi.ProcProcConeTolEnergy     = 'Tolerance energy for post-processing of cone';

%% layout of the panel
nrow        = 4;
rspacing    = 0.03;
ncol        = 3;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent handle of HEIDI panel children

h.qsm.panel.HEIDI = uipanel(hParent,...
        'Title','Homogeneity Enabled Incremental Dipole Inversion (HEIDI)',...
        'position',position,...
        'backgroundcolor',get(h.fig,'color'),'Visible','off');

%% Children of HEIDI panel
    
    panelParent = h.qsm.panel.HEIDI;

    % width of each element in a functional column, in normalised unit
    wratio = 0.65;
    
    % col 1, row 1
    % text|edit field pair: tolerance
    [h.qsm.HEIDI.text.tol,h.qsm.HEIDI.edit.tol] = sepia_construct_text_edit(...
        panelParent,'Tolerance:',               defaultTol,         [left(1) bottom(1) width height], wratio);
    
    % col 1, row 2
    % text|edit field pair: maximum iterations
    [h.qsm.HEIDI.text.maxIter,h.qsm.HEIDI.edit.maxIter] = sepia_construct_text_edit(...
        panelParent,'Max. iterations:',         defaultMaxIter,     [left(1) bottom(2) width height], wratio);
    
    % col 1, row 3
    % text|popup field pair: offset
    h.qsm.HEIDI.checkbox.isOffsetUse = uicontrol('Parent',h.qsm.panel.HEIDI ,...
        'Style','checkbox','String','Automatic offset determination',...
        'units','normalized','Position',[left(1) bottom(3) width height],...
        'backgroundcolor',get(h.fig,'color'),'Value',defaultOffsetUseBool);

    % col 1, row 4
    % text|popup field pair: offset
    h.qsm.HEIDI.checkbox.isFourierDomainFormula = uicontrol('Parent',h.qsm.panel.HEIDI ,...
        'Style','checkbox','String','Fourier domain formula',...
        'units','normalized','Position',[left(1) bottom(4) width height],...
        'backgroundcolor',get(h.fig,'color'),'Value',defaultIsFourierDomainFormula);

    % col 2, row 1
    % text|popup field pair: Tikhonov
    [h.qsm.HEIDI.text.Tikhonov,h.qsm.HEIDI.popup.Tikhonov] = sepia_construct_text_popup(...
        panelParent,'Tikhonov:',                menuTikhonov,	[left(2) bottom(1) width height], wratio);

    % col 2, row 2
    % text|popup field pair: Solver
    [h.qsm.HEIDI.text.solver,h.qsm.HEIDI.popup.solver] = sepia_construct_text_popup(...
        panelParent,'Solver:',                  menuSolver,	        [left(2) bottom(2) width height], wratio);

    % col 2, row 3
     % text|popup field pair: Solver
    [h.qsm.HEIDI.text.DipoleFilter,h.qsm.HEIDI.popup.DipoleFilter] = sepia_construct_text_popup(...
        panelParent,'Dipole Filter:',           menuDipoleFilter,   [left(2) bottom(3) width height], wratio);

    % col 3, row 1
    % text|edit field pair: Residual weighting
    [h.qsm.HEIDI.text.ResidualWeighting,h.qsm.HEIDI.edit.ResidualWeighting] = sepia_construct_text_edit(...
        panelParent,'Residual weighting:',     defaultResidualWeighting,        [left(3) bottom(1) width height], wratio);
    
    % col 3, row 2
    % text|edit field pair: Post-processing of Cone theshold
    [h.qsm.HEIDI.text.ProcProcConeThreshold,h.qsm.HEIDI.edit.ProcProcConeThreshold] = sepia_construct_text_edit(...
        panelParent,'Cone post-processing theshold:',     defaultProcProcConeThreshold,        [left(3) bottom(2) width height], wratio);

    % col 3, row 3
    % text|edit field pair: Post-processing of Cone tol
    [h.qsm.HEIDI.text.ProcProcConeTol,h.qsm.HEIDI.edit.ProcProcConeTol] = sepia_construct_text_edit(...
        panelParent,'Cone post-processing tolerance:',     defaultProcProcConeTol,        [left(3) bottom(3) width height], wratio);

    % col 3, row 4
    % text|edit field pair: Post-processing of Cone tol
    [h.qsm.HEIDI.text.ProcProcConeTolEnergy,h.qsm.HEIDI.edit.ProcProcConeTolEnergy] = sepia_construct_text_edit(...
        panelParent,'Cone post-processing energy:',     defaultProcProcConeTolEnergy,        [left(3) bottom(4) width height], wratio);
    
    
%% set tooltips
set(h.qsm.HEIDI.text.tol,                       'Tooltip',tooltip.qsm.heidi.tol);
set(h.qsm.HEIDI.text.maxIter,                   'Tooltip',tooltip.qsm.heidi.MaxIter);
set(h.qsm.HEIDI.checkbox.isOffsetUse,           'Tooltip',tooltip.qsm.heidi.isOffsetUse);
set(h.qsm.HEIDI.checkbox.isFourierDomainFormula,'Tooltip',tooltip.qsm.heidi.isFourierDomainFormula);
set(h.qsm.HEIDI.text.DipoleFilter,              'Tooltip',tooltip.qsm.heidi.DipoleFilter);
set(h.qsm.HEIDI.text.Tikhonov,                  'Tooltip',tooltip.qsm.heidi.Tikhonov);
set(h.qsm.HEIDI.text.solver,                    'Tooltip',tooltip.qsm.heidi.solver);
set(h.qsm.HEIDI.text.ResidualWeighting,         'Tooltip',tooltip.qsm.heidi.ResidualWeighting);
set(h.qsm.HEIDI.text.ProcProcConeThreshold,     'Tooltip',tooltip.qsm.heidi.ProcProcConeThreshold);
set(h.qsm.HEIDI.text.ProcProcConeTol,           'Tooltip',tooltip.qsm.heidi.ProcProcConeTol);
set(h.qsm.HEIDI.text.ProcProcConeTolEnergy,     'Tooltip',tooltip.qsm.heidi.ProcProcConeTolEnergy);

%% set callbacks
set(h.qsm.HEIDI.edit.maxIter,	                'Callback', {@EditInputMinMax_Callback,defaultMaxIter,  1,1});
set(h.qsm.HEIDI.edit.ResidualWeighting,     	'Callback', {@EditInputMinMax_Callback,defaultResidualWeighting,       0,0});
set(h.qsm.HEIDI.edit.ProcProcConeThreshold,     'Callback', {@EditInputMinMax_Callback,defaultProcProcConeThreshold,      0,0});
set(h.qsm.HEIDI.edit.ProcProcConeTol,    	    'Callback', {@EditInputMinMax_Callback,defaultProcProcConeTol,      0,0});
set(h.qsm.HEIDI.edit.ProcProcConeTolEnergy,     'Callback', {@EditInputMinMax_Callback,defaultProcProcConeTolEnergy,     0,0});

end
