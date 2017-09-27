% function qsm_hub
% gui
clear 

qsm_hub_AddPath

global h

fig=figure('Units','pixels','position',[400 20 500 650],...
    'MenuBar','None','Toolbar','None','Name','QSM hub','NumberTitle','off');

%% I/O
h.panel_dataIO = uipanel(fig,'Title','I/O',...
    'Position',[0.01 0.8 0.95 0.2],...
    'backgroundcolor',get(fig,'color'));
    h.text_input = uicontrol(h.panel_dataIO,'Style','text','String','Input dir:',...
        'Units','normalized','Position', [0.01 0.8 0.15 0.15],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(fig,'color'),...
        'tooltip','Input DICOM directory');
    h.edit_input = uicontrol('Parent',h.panel_dataIO,'Style','edit',...
        'units','normalized','position',[0.16 0.8 0.7 0.15],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    h.button_input = uicontrol('Parent',h.panel_dataIO,'Style','pushbutton','String','open',...
        'units','normalized','position',[0.87 0.8 0.1 0.15],...
        'backgroundcolor','white');
    h.text_output = uicontrol('Parent',h.panel_dataIO,'Style','text','String','Output dir:',...
        'units','normalized','Position',[0.01 0.64 0.15 0.15],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(fig,'color'),...
        'tooltip','Output directory');
    h.edit_output = uicontrol('Parent',h.panel_dataIO,'Style','edit',...
        'units','normalized','position',[0.16 0.64 0.7 0.15],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    h.button_output = uicontrol('Parent',h.panel_dataIO,'Style','pushbutton','String','open',...
        'units','normalized','position',[0.87 0.64 0.1 0.15],...
        'backgroundcolor','white');
    h.checkbox_brainExtraction = uicontrol('Parent',h.panel_dataIO,'Style','checkbox','String','FSL brain extraction',...
        'units','normalized','Position',[0.01 0.4 1 0.2],...
        'backgroundcolor',get(fig,'color'),...
        'tooltip','Use FSL brain extraction (bet)');
    h.text_maskdir = uicontrol('Parent',h.panel_dataIO,'Style','text','String','Mask dir:',...
        'units','normalized','Position',[0.01 0.24 0.15 0.15],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(fig,'color'),...
        'tooltip','Specify mask file (NIfTI_GZ format)');
    h.edit_maskdir = uicontrol('Parent',h.panel_dataIO,'Style','edit',...
        'units','normalized','position',[0.16 0.24 0.7 0.15],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    h.button_maskdir = uicontrol('Parent',h.panel_dataIO,'Style','pushbutton','String','open',...
        'units','normalized','position',[0.87 0.24 0.1 0.15],...
        'backgroundcolor','white');

%% phase unwrap
h.panel_phaseUnwrap = uipanel(fig,'Title','Total field recovery and phase unwrapping',...
    'position',[0.01 0.59 0.95 0.2]);
    h.text_phaseUnwrap = uicontrol('Parent',h.panel_phaseUnwrap,'Style','text','String','Phase unwrapping:',...
        'units','normalized','position',[0.01 0.84 0.3 0.15],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(fig,'color'),...
        'tooltip','Select phase unwrapping algorithm');
    h.popup_phaseUnwrap = uicontrol('Parent',h.panel_phaseUnwrap,'Style','popup',...
        'String',{'Laplacian','Jena','Region growing','Graphcut'},...
        'units','normalized','position',[0.31 0.84 0.4 0.15]); 
    h.text_unit = uicontrol('Parent',h.panel_phaseUnwrap,'Style','text','String','Output unit:',...
        'units','normalized','position',[0.01 0.65 0.3 0.15],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(fig,'color'),...
        'tooltip','Select unwrapped phase map unit');
    h.popup_unit = uicontrol('Parent',h.panel_phaseUnwrap,'Style','popup',...
        'String',{'radHz','ppm','rad','Hz'},...
        'units','normalized','position',[0.31 0.65 0.4 0.15]);
    
%% background field
h.panel_bkgRemoval = uipanel(fig,'Title','Background field removal',...
    'position',[0.01 0.33 0.95 0.25],...
    'backgroundcolor',get(fig,'color'));
    h.text_bkgRemoval = uicontrol('Parent',h.panel_bkgRemoval,'Style','text','String','Method:',...
        'units','normalized','position',[0.01 0.85 0.3 0.1],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(fig,'color'),...
        'tooltip','Select background field removal method');
    h.popup_bkgRemoval = uicontrol('Parent',h.panel_bkgRemoval,'Style','popup',...
        'String',{'LBV','PDF','RESHARP','SHARP','VSHARP STI suite','VSHARP','iHARPERELLA'},...
        'units','normalized','position',[0.31 0.85 0.4 0.1]) ;   
    h.checkbox_bkgRemoval = uicontrol('Parent',h.panel_bkgRemoval,'Style','checkbox','String','Refine local field by polynomial fit',...
        'units','normalized','position',[0.01 0.01 1 0.1],...
        'backgroundcolor',get(fig,'color'));
    
    % LBV
    h.panel_bkgRemoval_LBV = uipanel(h.panel_bkgRemoval,'Title','Laplacian boundary value (LBV)',...
        'position',[0.01 0.15 0.95 0.65],...
        'backgroundcolor',get(fig,'color'));
        h.text_LBV_tol = uicontrol('Parent',h.panel_bkgRemoval_LBV,'Style','text',...
            'String','Tolerance:',...
            'units','normalized','position',[0.01 0.75 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'),...
            'tooltip','Iteration stopping criteria on te coarsest grid');
        h.edit_LBV_tol = uicontrol('Parent',h.panel_bkgRemoval_LBV,'Style','edit',...
            'String','0.0001',...
            'units','normalized','position',[0.25 0.75 0.2 0.2],...
            'backgroundcolor','white');
        h.text_LBV_depth = uicontrol('Parent',h.panel_bkgRemoval_LBV,'Style','text',...
            'String','Depth:',...
            'units','normalized','position',[0.01 0.5 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'),...
            'tooltip','Multigrid level; The largest length scale is 2^depth * voxel size');
        h.edit_LBV_depth = uicontrol('Parent',h.panel_bkgRemoval_LBV,'Style','edit',...
            'String','2',...
            'units','normalized','position',[0.25 0.5 0.2 0.2],...
            'backgroundcolor','white');
        h.text_LBV_peel = uicontrol('Parent',h.panel_bkgRemoval_LBV,'Style','text',...
            'String','Peel:',...
            'units','normalized','position',[0.01 0.25 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'),...
            'tooltip','No. of boundary layers to be peeled off');
        h.edit_LBV_peel = uicontrol('Parent',h.panel_bkgRemoval_LBV,'Style','edit',...
            'String','5',...
            'units','normalized','position',[0.25 0.25 0.2 0.2],...
            'backgroundcolor','white');
        
    % PDF
    h.panel_bkgRemoval_PDF = uipanel(h.panel_bkgRemoval,'Title','Projection onto dipole field (PDF)',...
        'position',[0.01 0.15 0.95 0.65],...
        'backgroundcolor',get(fig,'color'),'Visible','off');
        h.text_PDF_tol = uicontrol('Parent',h.panel_bkgRemoval_PDF,'Style','text',...
            'String','Tolerance:',...
            'units','normalized','position',[0.01 0.75 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'));
        h.edit_PDF_tol = uicontrol('Parent',h.panel_bkgRemoval_PDF,'Style','edit',...
            'String','0.0001',...
            'units','normalized','position',[0.25 0.75 0.2 0.2],...
            'backgroundcolor','white');
        h.text_PDF_maxIter = uicontrol('Parent',h.panel_bkgRemoval_PDF,'Style','text',...
            'String','Iterations:',...
            'units','normalized','position',[0.01 0.5 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'),...
            'tooltip','Maximum iterations');
        h.edit_PDF_maxIter = uicontrol('Parent',h.panel_bkgRemoval_PDF,'Style','edit',...
            'String','50',...
            'units','normalized','position',[0.25 0.5 0.2 0.2],...
            'backgroundcolor','white');
        h.text_PDF_cgSolver = uicontrol('Parent',h.panel_bkgRemoval_PDF,'Style','text',...
            'String','CG solver',...
            'units','normalized','position',[0.01 0.25 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'),...
            'tooltip','Select CG solver');
        h.popup_PDF_cgSolver = uicontrol('Parent',h.panel_bkgRemoval_PDF,'Style','popup',...
            'String',{'MEDI cgsolver','Matlab pcg'},...
            'units','normalized','position',[0.25 0.25 0.5 0.2]);
        
    % SHARP
    h.panel_bkgRemoval_SHARP = uipanel(h.panel_bkgRemoval,'Title','SHARP',...
        'position',[0.01 0.15 0.95 0.65],...
        'backgroundcolor',get(fig,'color'),'Visible','off');
        h.text_SHARP_radius = uicontrol('Parent',h.panel_bkgRemoval_SHARP,'Style','text',...
            'String','Radius:',...
            'units','normalized','position',[0.01 0.75 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'),...
            'tooltip','Radius of sphere; unit in voxel');
        h.edit_SHARP_radius = uicontrol('Parent',h.panel_bkgRemoval_SHARP,'Style','edit',...
            'String','4',...
            'units','normalized','position',[0.25 0.75 0.2 0.2],...
            'backgroundcolor','white');
        h.text_SHARP_threshold = uicontrol('Parent',h.panel_bkgRemoval_SHARP,'Style','text',...
            'String','Threshold:',...
            'units','normalized','position',[0.01 0.5 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'),...
            'tooltip','Threshold used in truncated SVD');
        h.edit_SHARP_threshold = uicontrol('Parent',h.panel_bkgRemoval_SHARP,'Style','edit',...
            'String','0.03',...
            'units','normalized','position',[0.25 0.5 0.2 0.2],...
            'backgroundcolor','white');
    
    % RESHARP    
    h.panel_bkgRemoval_RESHARP = uipanel(h.panel_bkgRemoval,'Title','Regularisation SHARP',...
        'position',[0.01 0.15 0.95 0.65],...
        'backgroundcolor',get(fig,'color'),'Visible','off');
        h.text_RESHARP_radius = uicontrol('Parent',h.panel_bkgRemoval_RESHARP,'Style','text',...
            'String','Radius:',...
            'units','normalized','position',[0.01 0.75 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'));
        h.edit_RESHARP_radius = uicontrol('Parent',h.panel_bkgRemoval_RESHARP,'Style','edit',...
            'String','4',...
            'units','normalized','position',[0.25 0.75 0.2 0.2],...
            'backgroundcolor','white');
        h.text_RESHARP_lambda = uicontrol('Parent',h.panel_bkgRemoval_RESHARP,'Style','text',...
            'String','Regularisation:',...
            'units','normalized','position',[0.01 0.5 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'));
        h.edit_RESHARP_lambda = uicontrol('Parent',h.panel_bkgRemoval_RESHARP,'Style','edit',...
            'String','0.01',...
            'units','normalized','position',[0.25 0.5 0.2 0.2],...
            'backgroundcolor','white');
        
    % VSHARPSTI    
    h.panel_bkgRemoval_VSHARPSTI = uipanel(h.panel_bkgRemoval,'Title','STI suite Variable SHARP',...
        'position',[0.01 0.15 0.95 0.65],...
        'backgroundcolor',get(fig,'color'),'Visible','off');
    
    % VSHARP
    h.panel_bkgRemoval_VSHARP = uipanel(h.panel_bkgRemoval,'Title','Variable SHARP',...
        'position',[0.01 0.15 0.95 0.65],...
        'backgroundcolor',get(fig,'color'),'Visible','off');
        h.text_VSHARP_maxRadius = uicontrol('Parent',h.panel_bkgRemoval_VSHARP,'Style','text',...
            'String','Max. radius:',...
            'units','normalized','position',[0.01 0.75 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'));
        h.edit_VSHARP_maxRadius = uicontrol('Parent',h.panel_bkgRemoval_VSHARP,'Style','edit',...
            'String','10',...
            'units','normalized','position',[0.25 0.75 0.2 0.2],...
            'backgroundcolor','white');
        h.text_VSHARP_minRadius = uicontrol('Parent',h.panel_bkgRemoval_VSHARP,'Style','text',...
            'String','Min. radius:',...
            'units','normalized','position',[0.01 0.5 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'));
        h.edit_VSHARP_minRadius = uicontrol('Parent',h.panel_bkgRemoval_VSHARP,'Style','edit',...
            'String','3',...
            'units','normalized','position',[0.25 0.5 0.2 0.2],...
            'backgroundcolor','white');
        
    % iHARPERELLA    
    h.panel_bkgRemoval_iHARPERELLA = uipanel(h.panel_bkgRemoval,'Title','iHARPERELLA',...
        'position',[0.01 0.15 0.95 0.65],...
        'backgroundcolor',get(fig,'color'),'Visible','off');
        h.text_iHARPERELLA_maxIter = uicontrol('Parent',h.panel_bkgRemoval_iHARPERELLA,'Style','text',...
            'String','Max. iterations:',...
            'units','normalized','position',[0.01 0.75 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'));
        h.edit_iHARPERELLA_maxIter = uicontrol('Parent',h.panel_bkgRemoval_iHARPERELLA,'Style','edit',...
            'String','100',...
            'units','normalized','position',[0.25 0.75 0.2 0.2],...
            'backgroundcolor','white');

%% QSM
h.panel_qsm = uipanel(fig,'Title','QSM','backgroundcolor',get(fig,'color'),...
    'position',[0.01 0.07 0.95 0.25]);
    h.text_qsm = uicontrol('Parent',h.panel_qsm,'Style','text','String','Method:',...
        'units','normalized','position',[0.01 0.85 0.15 0.1],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(fig,'color'));
    h.popup_qsm = uicontrol('Parent',h.panel_qsm,'Style','popup',...
        'String',{'TKD','Closed-form solution','STI suite iLSQR','iLSQR','FANSI'},...
        'units','normalized','position',[0.31 0.85 0.4 0.1]) ; 
    
    % TKD    
    h.panel_qsm_TKD = uipanel(h.panel_qsm,'Title','Thresholded k-space division (TKD)',...
        'position',[0.01 0.04 0.95 0.75],...
        'backgroundcolor',get(fig,'color'),'Visible','on');
        h.text_TKD_threshold = uicontrol('Parent',h.panel_qsm_TKD,'Style','text',...
            'String','Threshold:',...
            'units','normalized','position',[0.01 0.75 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'));
        h.edit_TKD_threshold = uicontrol('Parent',h.panel_qsm_TKD,'Style','edit',...
            'String','0.15',...
            'units','normalized','position',[0.25 0.75 0.2 0.2],...
            'backgroundcolor','white');

    % Closed-form solution  
    h.panel_qsm_cfs = uipanel(h.panel_qsm,'Title','Closed-form solution with L2-norm regularisation',...
        'position',[0.01 0.04 0.95 0.75],...
        'backgroundcolor',get(fig,'color'),'Visible','off');
        h.text_cfs_lambda = uicontrol('Parent',h.panel_qsm_cfs,'Style','text',...
            'String','Lambda:',...
            'units','normalized','position',[0.01 0.75 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'));
        h.edit_cfs_lambda = uicontrol('Parent',h.panel_qsm_cfs,'Style','edit',...
            'String','0.13',...
            'units','normalized','position',[0.25 0.75 0.2 0.2],...
            'backgroundcolor','white');
        h.checkbox_cfs_lambda = uicontrol('Parent',h.panel_qsm_cfs,'Style','checkbox',...
            'String','Self-optimisation by L-curve approach',...
            'units','normalized','position',[0.01 0.5 1 0.2],...
            'backgroundcolor',get(fig,'color'));
    
    % iLSQR  
    h.panel_qsm_iLSQR = uipanel(h.panel_qsm,'Title','Iterative LSQR',...
        'position',[0.01 0.04 0.95 0.75],...
        'backgroundcolor',get(fig,'color'),'Visible','off');
        h.text_iLSQR_tol = uicontrol('Parent',h.panel_qsm_iLSQR,'Style','text',...
            'String','Tolerance:',...
            'units','normalized','position',[0.01 0.75 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'));
        h.edit_iLSQR_tol = uicontrol('Parent',h.panel_qsm_iLSQR,'Style','edit',...
            'String','0.001',...
            'units','normalized','position',[0.25 0.75 0.2 0.2],...
            'backgroundcolor','white');
        h.text_iLSQR_maxIter = uicontrol('Parent',h.panel_qsm_iLSQR,'Style','text',...
            'String','Iterations:',...
            'units','normalized','position',[0.01 0.5 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'));
        h.edit_iLSQR_maxIter = uicontrol('Parent',h.panel_qsm_iLSQR,'Style','edit',...
            'String','100',...
            'units','normalized','position',[0.25 0.5 0.2 0.2],...
            'backgroundcolor','white');
        h.text_iLSQR_lambda = uicontrol('Parent',h.panel_qsm_iLSQR,'Style','text',...
            'String','Lambda:',...
            'units','normalized','position',[0.01 0.25 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'));
        h.edit_iLSQR_lambda = uicontrol('Parent',h.panel_qsm_iLSQR,'Style','edit',...
            'String','0.13',...
            'units','normalized','position',[0.25 0.25 0.2 0.2],...
            'backgroundcolor','white');
        h.checkbox_iLSQR_lambda = uicontrol('Parent',h.panel_qsm_iLSQR,'Style','checkbox',...
            'String','Self-optimisation by L-curve approach',...
            'units','normalized','position',[0.01 0.01 1 0.2],...
            'backgroundcolor',get(fig,'color'));
        
    % STI suite iLSQR  
    h.panel_qsm_STIiLSQR = uipanel(h.panel_qsm,'Title','STI suite iLSQR',...
        'position',[0.01 0.04 0.95 0.75],...
        'backgroundcolor',get(fig,'color'),'Visible','off');
        h.text_STIiLSQR_threshold = uicontrol('Parent',h.panel_qsm_STIiLSQR,'Style','text',...
            'String','Threshold:',...
            'units','normalized','position',[0.01 0.75 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'));
        h.edit_STIiLSQR_threshold = uicontrol('Parent',h.panel_qsm_STIiLSQR,'Style','edit',...
            'String','0.01',...
            'units','normalized','position',[0.25 0.75 0.2 0.2],...
            'backgroundcolor','white');
        h.text_STIiLSQR_maxIter = uicontrol('Parent',h.panel_qsm_STIiLSQR,'Style','text',...
            'String','Iterations:',...
            'units','normalized','position',[0.01 0.5 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'));
        h.edit_STIiLSQR_maxIter = uicontrol('Parent',h.panel_qsm_STIiLSQR,'Style','edit',...
            'String','100',...
            'units','normalized','position',[0.25 0.5 0.2 0.2],...
            'backgroundcolor','white');
        h.text_STIiLSQR_tol1 = uicontrol('Parent',h.panel_qsm_STIiLSQR,'Style','text',...
            'String','Tolerance 1:',...
            'units','normalized','position',[0.01 0.25 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'));
        h.edit_STIiLSQR_tol1 = uicontrol('Parent',h.panel_qsm_STIiLSQR,'Style','edit',...
            'String','0.01',...
            'units','normalized','position',[0.25 0.25 0.2 0.2],...
            'backgroundcolor','white');
        h.text_STIiLSQR_tol2 = uicontrol('Parent',h.panel_qsm_STIiLSQR,'Style','text',...
            'String','Tolerance 2:',...
            'units','normalized','position',[0.5 0.25 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'));
        h.edit_STIiLSQR_tol2 = uicontrol('Parent',h.panel_qsm_STIiLSQR,'Style','edit',...
            'String','0.001',...
            'units','normalized','position',[0.75 0.25 0.2 0.2],...
            'backgroundcolor','white');
        h.text_STIiLSQR_padSize = uicontrol('Parent',h.panel_qsm_STIiLSQR,'Style','text',...
            'String','Pad size:',...
            'units','normalized','position',[0.01 0.01 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'));
        h.edit_STIiLSQR_padSize = uicontrol('Parent',h.panel_qsm_STIiLSQR,'Style','edit',...
            'String','4',...
            'units','normalized','position',[0.25 0.01 0.2 0.2],...
            'backgroundcolor','white');
    
    % FANSI
    h.panel_qsm_FANSI = uipanel(h.panel_qsm,'Title','FANSI',...
        'position',[0.01 0.04 0.95 0.75],...
        'backgroundcolor',get(fig,'color'),'Visible','off');
        h.text_FANSI_tol = uicontrol('Parent',h.panel_qsm_FANSI,'Style','text',...
            'String','Tolerance:',...
            'units','normalized','position',[0.01 0.75 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'));
        h.edit_FANSI_tol = uicontrol('Parent',h.panel_qsm_FANSI,'Style','edit',...
            'String','1',...
            'units','normalized','position',[0.25 0.75 0.2 0.2],...
            'backgroundcolor','white');
        h.text_FANSI_lambda = uicontrol('Parent',h.panel_qsm_FANSI,'Style','text',...
            'String','L1 penalty:',...
            'units','normalized','position',[0.01 0.5 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'));
        h.edit_FANSI_lambda = uicontrol('Parent',h.panel_qsm_FANSI,'Style','edit',...
            'String','3e-5',...
            'units','normalized','position',[0.25 0.5 0.2 0.2],...
            'backgroundcolor','white');
        h.text_FANSI_mu = uicontrol('Parent',h.panel_qsm_FANSI,'Style','text',...
            'String','Consistency:',...
            'units','normalized','position',[0.5 0.5 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'));
        h.edit_FANSI_mu = uicontrol('Parent',h.panel_qsm_FANSI,'Style','edit',...
            'String','5e-5',...
            'units','normalized','position',[0.75 0.5 0.2 0.2],...
            'backgroundcolor','white');
        h.text_FANSI_maxIter = uicontrol('Parent',h.panel_qsm_FANSI,'Style','text',...
            'String','Max. iterations:',...
            'units','normalized','position',[0.01 0.25 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'));
        h.edit_FANSI_maxIter = uicontrol('Parent',h.panel_qsm_FANSI,'Style','edit',...
            'String','50',...
            'units','normalized','position',[0.25 0.25 0.2 0.2],...
            'backgroundcolor','white');
        h.text_FANSI_solver = uicontrol('Parent',h.panel_qsm_FANSI,'Style','text',...
            'String','Solver:',...
            'units','normalized','position',[0.01 0.01 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'));
        h.popup_FANSI_solver = uicontrol('Parent',h.panel_qsm_FANSI,'Style','popup',...
            'String',{'Linear','Non-linear'},...
            'units','normalized','position',[0.25 0.01 0.2 0.2],...
            'backgroundcolor','white');
        h.text_FANSI_constraints = uicontrol('Parent',h.panel_qsm_FANSI,'Style','text',...
            'String','Constraint:',...
            'units','normalized','position',[0.5 0.01 0.2 0.2],...
            'HorizontalAlignment','left',...
            'backgroundcolor',get(fig,'color'));
        h.popup_FANSI_constraints = uicontrol('Parent',h.panel_qsm_FANSI,'Style','popup',...
            'String',{'TV','TGV'},...
            'units','normalized','position',[0.75 0.01 0.2 0.2],...
            'backgroundcolor','white');

h.pushbutton_start = uicontrol('Parent',fig,'Style','pushbutton',...
    'String','Start',...
    'units','normalized','Position',[0.85 0.01 0.1 0.05],...
    'backgroundcolor',get(fig,'color'));
    
%% Set Callback
set(h.button_input,           	'Callback',{@ButtonGetInputDir_Callback});
set(h.button_output,            'Callback',{@ButtonGetOutputDir_Callback});
set(h.checkbox_brainExtraction,	'Callback',{@CheckboxBrainExtraction_Callback});
set(h.button_maskdir,       	'Callback',{@ButtonGetMaskDir_Callback});
% set(h.popup_phaseUnwrap,          'Callback',{@imageselect_Callback,h});
% set(h.popup_unit,                 'Callback',{@imageselect_Callback,h});
set(h.popup_bkgRemoval,      	'Callback',{@PopupBkgRemoval_Callback});
set(h.popup_qsm,                'Callback',{@PopupQSM_Callback});
set(h.checkbox_cfs_lambda,      'Callback',{@CheckboxCFS_Callback});
set(h.checkbox_iLSQR_lambda,    'Callback',{@CheckboxiLSQR_Callback});
set(h.pushbutton_start,           'Callback',{@PushbuttonStart_Callback});
% end
%% Callback
function ButtonGetInputDir_Callback(source,eventdata)

global h

pathDir = uigetdir;

if pathDir ~= 0
    set(h.edit_input,'String',pathDir);
    parts = strfind(pathDir, '/');
    pathDirParent = pathDir(1:parts(end));
    set(h.edit_output,'String',[pathDirParent 'output']);
end
end

function ButtonGetOutputDir_Callback(source,eventdata)

global h

pathDir = uigetdir;

if pathDir ~= 0
    set(h.edit_output,'String',pathDir);
end
end

function CheckboxBrainExtraction_Callback(source,eventdata)

global h

if ~h.checkbox_brainExtraction.Value
    set(h.button_maskdir,'Enable','on');
    set(h.edit_maskdir,'Enable','on');
else
    set(h.button_maskdir,'Enable','off');
    set(h.edit_maskdir,'Enable','off');
end
end
    
function ButtonGetMaskDir_Callback(source,eventdata)

global h

pathDir = uigetfile({'*.nii.gz';'*.nii'},'Select the mask file');

if pathDir ~= 0
    set(h.edit_maskdir,'String',pathDir);
end
end

function PopupBkgRemoval_Callback(source,eventdata)

global h

method = source.String{source.Value,1} ;

switch method
    case 'LBV'
        set(h.panel_bkgRemoval_LBV,'Visible','on');
        set(h.panel_bkgRemoval_PDF,'Visible','off');
        set(h.panel_bkgRemoval_RESHARP,'Visible','off');
        set(h.panel_bkgRemoval_SHARP,'Visible','off');
        set(h.panel_bkgRemoval_VSHARP,'Visible','off');
        set(h.panel_bkgRemoval_VSHARPSTI,'Visible','off');
        set(h.panel_bkgRemoval_iHARPERELLA,'Visible','off');
    case 'PDF'
        set(h.panel_bkgRemoval_LBV,'Visible','off');
        set(h.panel_bkgRemoval_PDF,'Visible','on');
        set(h.panel_bkgRemoval_RESHARP,'Visible','off');
        set(h.panel_bkgRemoval_SHARP,'Visible','off');
        set(h.panel_bkgRemoval_VSHARP,'Visible','off');
        set(h.panel_bkgRemoval_VSHARPSTI,'Visible','off');
        set(h.panel_bkgRemoval_iHARPERELLA,'Visible','off');
    case 'RESHARP'
        set(h.panel_bkgRemoval_LBV,'Visible','off');
        set(h.panel_bkgRemoval_PDF,'Visible','off');
        set(h.panel_bkgRemoval_RESHARP,'Visible','on');
        set(h.panel_bkgRemoval_SHARP,'Visible','off');
        set(h.panel_bkgRemoval_VSHARP,'Visible','off');
        set(h.panel_bkgRemoval_VSHARPSTI,'Visible','off');
        set(h.panel_bkgRemoval_iHARPERELLA,'Visible','off');
    case 'SHARP'
        set(h.panel_bkgRemoval_LBV,'Visible','off');
        set(h.panel_bkgRemoval_PDF,'Visible','off');
        set(h.panel_bkgRemoval_RESHARP,'Visible','off');
        set(h.panel_bkgRemoval_SHARP,'Visible','on');
        set(h.panel_bkgRemoval_VSHARP,'Visible','off');
        set(h.panel_bkgRemoval_VSHARPSTI,'Visible','off');
        set(h.panel_bkgRemoval_iHARPERELLA,'Visible','off');
    case 'VSHARP'
        set(h.panel_bkgRemoval_LBV,'Visible','off');
        set(h.panel_bkgRemoval_PDF,'Visible','off');
        set(h.panel_bkgRemoval_RESHARP,'Visible','off');
        set(h.panel_bkgRemoval_SHARP,'Visible','off');
        set(h.panel_bkgRemoval_VSHARP,'Visible','on');
        set(h.panel_bkgRemoval_VSHARPSTI,'Visible','off');
        set(h.panel_bkgRemoval_iHARPERELLA,'Visible','off');
    case 'VSHARP STI suite'
        set(h.panel_bkgRemoval_LBV,'Visible','off');
        set(h.panel_bkgRemoval_PDF,'Visible','off');
        set(h.panel_bkgRemoval_RESHARP,'Visible','off');
        set(h.panel_bkgRemoval_SHARP,'Visible','off');
        set(h.panel_bkgRemoval_VSHARP,'Visible','off');
        set(h.panel_bkgRemoval_VSHARPSTI,'Visible','on');
        set(h.panel_bkgRemoval_iHARPERELLA,'Visible','off');
    case 'iHARPERELLA'
        set(h.panel_bkgRemoval_LBV,'Visible','off');
        set(h.panel_bkgRemoval_PDF,'Visible','off');
        set(h.panel_bkgRemoval_RESHARP,'Visible','off');
        set(h.panel_bkgRemoval_SHARP,'Visible','off');
        set(h.panel_bkgRemoval_VSHARP,'Visible','off');
        set(h.panel_bkgRemoval_VSHARPSTI,'Visible','off');
        set(h.panel_bkgRemoval_iHARPERELLA,'Visible','on');
end
end

function PopupQSM_Callback(source,eventdata)

global h

method = source.String{source.Value,1} ;

switch method
    case 'TKD'
        set(h.panel_qsm_TKD,'Visible','on');
        set(h.panel_qsm_cfs,'Visible','off');
        set(h.panel_qsm_iLSQR,'Visible','off');
        set(h.panel_qsm_STIiLSQR,'Visible','off');
        set(h.panel_qsm_FANSI,'Visible','off');
    case 'Closed-form solution'
        set(h.panel_qsm_TKD,'Visible','off');
        set(h.panel_qsm_cfs,'Visible','on');
        set(h.panel_qsm_iLSQR,'Visible','off');
        set(h.panel_qsm_STIiLSQR,'Visible','off');
        set(h.panel_qsm_FANSI,'Visible','off');
    case 'STI suite iLSQR'
        set(h.panel_qsm_TKD,'Visible','off');
        set(h.panel_qsm_cfs,'Visible','off');
        set(h.panel_qsm_iLSQR,'Visible','off');
        set(h.panel_qsm_STIiLSQR,'Visible','on');
        set(h.panel_qsm_FANSI,'Visible','off');
    case 'iLSQR'
        set(h.panel_qsm_TKD,'Visible','off');
        set(h.panel_qsm_cfs,'Visible','off');
        set(h.panel_qsm_iLSQR,'Visible','on');
        set(h.panel_qsm_STIiLSQR,'Visible','off');
        set(h.panel_qsm_FANSI,'Visible','off');
    case 'FANSI'
        set(h.panel_qsm_TKD,'Visible','off');
        set(h.panel_qsm_cfs,'Visible','off');
        set(h.panel_qsm_iLSQR,'Visible','off');
        set(h.panel_qsm_STIiLSQR,'Visible','off');
        set(h.panel_qsm_FANSI,'Visible','on');
end
end

function CheckboxCFS_Callback(source,eventdata)

global h

if ~h.checkbox_cfs_lambda.Value
    set(h.edit_cfs_lambda,'Enable','on');
else
    set(h.edit_cfs_lambda,'Enable','off');
end
end

function CheckboxiLSQR_Callback(source,eventdata)

global h

if ~h.checkbox_iLSQR_lambda.Value
    set(h.edit_iLSQR_lambda,'Enable','on');
else
    set(h.edit_iLSQR_lambda,'Enable','off');
end
end

function PushbuttonStart_Callback(source,eventdata)

global h

subsampling=1;
BFR_tol=1e-4;BFR_depth=4;BFR_peel=2;BFR_iteration=50;
BFR_CGdefault=true;BFR_radius=4;BFR_alpha=0.01;BFR_threshold=0.03;
QSM_threshold=0.15;QSM_lambda=0.13;QSM_optimise=false;
QSM_tol=1e-3;QSM_maxiter=50;QSM_tol1=0.01;QSM_tol2=0.001;QSM_padsize=[4,4,4];
QSM_mu1=5e-5;QSM_solver='linear';QSM_constraint='tv';

inputDir = get(h.edit_input,'String');
outputDir = get(h.edit_output,'String');
maskDir = get(h.edit_maskdir,'String');
isBET = get(h.checkbox_brainExtraction,'Value');
phaseUnwrap = h.popup_phaseUnwrap.String{h.popup_phaseUnwrap.Value,1};
units = h.popup_unit.String{h.popup_unit.Value,1};
BFR = h.popup_bkgRemoval.String{h.popup_bkgRemoval.Value,1};
QSM_method = h.popup_qsm.String{h.popup_qsm.Value,1};
refine = get(h.checkbox_bkgRemoval,'Value');

% get backgroud field removal algorithm parameters
switch BFR
    case 'LBV'
        BFR='lbv';
        try BFR_tol = str2double(get(h.edit_LBV_tol,'String')); catch; BFR_tol=1e-4; end
        try BFR_depth = str2double(get(h.edit_LBV_depth,'String')); catch; BFR_depth=4; end
        try BFR_peel = str2double(get(h.edit_LBV_peel,'String')); catch; BFR_peel=4; end
    case 'PDF'
        BFR='pdf';
        try BFR_tol = str2double(get(h.edit_PDF_tol,'String')); catch; BFR_tol=1e-2; end
        try BFR_iteration = str2double(get(h.edit_PDF_maxIter,'String')); catch; BFR_iteration=50; end
        try BFR_CGdefault = h.popup_PDF_cgSolver.String{h.popup_PDF_cgSolver.Value,1}; catch; BFR_CGdefault=true; end
    case 'RESHARP'
        BFR='resharp';
        try BFR_radius = str2double(get(h.edit_RESHARP_radius,'String')); catch; BFR_radius=4; end
        try BFR_alpha = str2double(get(h.edit_RESHARP_lambda,'String')); catch; BFR_alpha=0.01; end
    case 'SHARP'
        BFR='sharp';
        try BFR_radius = str2double(get(h.edit_SHARP_radius,'String')); catch; BFR_radius=4; end
        try BFR_threshold = str2double(get(h.edit_SHARP_threshold,'String')); catch; BFR_threshold=0.03; end
    case 'VSHARP'
        BFR='vsharp';
        try maxRadius = str2double(get(h.edit_VSHARP_maxRadius,'String')); catch; maxRadius=10; end
        try minRadius = str2double(get(h.edit_VSHARP_minRadius,'String')); catch; minRadius=3; end
        BFR_radius = maxRadius:-2:minRadius;
    case 'iHARPERELLA'
        BFR='iharperella';
        try BFR_iteration = str2double(get(h.edit_iHARPERELLA_maxIter,'String')); catch; BFR_iteration=100; end 
    case 'VSHARP STI suite'
        BFR='vsharpsti';
end

% get QSM algorithm parameters
switch QSM_method
    case 'TKD'
        QSM_method='tkd';
        try QSM_threshold = str2double(get(h.edit_TKD_threshold,'String')); catch; QSM_threshold=0.15; end
    case 'Closed-form solution'
        QSM_method='closedforml2';
        try QSM_lambda = str2double(get(h.edit_cfs_lambda,'String')); catch; QSM_lambda=0.13; end
        try QSM_optimise = get(h.checkbox_cfs_lambda,'Value'); catch; QSM_optimise=false; end
    case 'STI suite iLSQR'
        QSM_method='stisuiteilsqr';
        try QSM_threshold = str2double(get(h.edit_STIiLSQR_threshold,'String')); catch; QSM_threshold=0.01; end
        try QSM_maxiter = str2double(get(h.edit_STIiLSQR_maxIter,'String')); catch; QSM_maxiter=100; end
        try QSM_tol1 = str2double(get(h.edit_STIiLSQR_tol1,'String')); catch; QSM_tol1=0.01; end
        try QSM_tol2 = str2double(get(h.edit_STIiLSQR_tol2,'String')); catch; QSM_tol2=0.001; end
        try QSM_padsize = str2double(get(h.edit_STIiLSQR_padSize,'String')); catch; QSM_padsize=4; end
        QSM_padsize = [QSM_padsize,QSM_padsize,QSM_padsize];
    case 'iLSQR'
        QSM_method='ilsqr';
        try QSM_tol = str2double(get(h.edit_iLSQR_tol,'String')); catch; QSM_tol=0.001; end
        try QSM_maxiter = str2double(get(h.edit_iLSQR_maxIter,'String')); catch; QSM_maxiter=100; end
        try QSM_lambda = str2double(get(h.edit_iLSQR_lambda,'String')); catch; QSM_lambda=0.13; end
        try QSM_optimise = get(h.checkbox_iLSQR_lambda,'Value'); catch; QSM_optimise=false; end 
    case 'FANSI'
        QSM_method='fansi';
        try QSM_tol = str2double(get(h.edit_FANSI_tol,'String')); catch; QSM_tol=1; end
        try QSM_lambda = str2double(get(h.edit_FANSI_lambda,'String')); catch; QSM_lambda=3e-5; end
        try QSM_mu1 = str2double(get(h.edit_FANSI_mu,'String')); catch; QSM_mu1=5e-5; end
        try QSM_maxiter = str2double(get(h.edit_FANSI_maxIter,'String')); catch; QSM_maxiter=50; end
        try 
            QSM_solver = h.popup_FANSI_solver.String{h.popup_FANSI_solver.Value,1}; 
        catch
            QSM_solver='linear'; 
        end 
        try 
            QSM_constraint = h.popup_FANSI_constraints.String{h.popup_FANSI_constraints.Value,1}; 
        catch
            QSM_constraint='tv'; 
        end 
end

% look for mask
try 
    mask = load_nii_img_only(maskDir);
catch
    mask = [];
end

QSMHub(inputDir,outputDir,'FSLBet',isBET,'mask',mask,'unwrap',phaseUnwrap,...
    'unit',units,'Subsampling',subsampling,'BFR',BFR,'refine',refine,'BFR_tol',BFR_tol,...
    'depth',BFR_depth,'peel',BFR_peel,'BFR_iteration',BFR_iteration,'CGsolver',BFR_CGdefault,...
    'BFR_radius',BFR_radius,'BFR_alpha',BFR_alpha,'BFR_threshold',BFR_threshold,...
    'QSM',QSM_method,'QSM_threshold',QSM_threshold,'QSM_lambda',QSM_lambda,'QSM_optimise',QSM_optimise,...
    'QSM_tol',QSM_tol,'QSM_iteration',QSM_maxiter,'QSM_tol1',QSM_tol1,'QSM_tol2',QSM_tol2,...
    'QSM_padsize',QSM_padsize,'QSM_mu',QSM_mu1,QSM_solver,QSM_constraint);

end

%% TODO: excludes unrealiable bkg remove