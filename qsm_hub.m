%% qsm_hub
%
% Description: This is a GUI of qsm_hub, which is a pipeline control tool
% for standard QSM processing. It supports the following processing steps:
% (1) phase unwrapping
% (2) background field removal
% (3) QSM
% All these processing steps can be performed in one go or worked as
% standalones.
%
% After perfoming all the selected processing steps sucessfully, a log file
% will be generated. This log file contains the text of how a function is
% called to run the process(es). You can also use the content of the log
% file to perform your data processing without calling this GUI  
% (which in some cases good for batch job).
%
% To run this GUI, simply add the directory containing this file in your
% Matlab's PATH. This GUI will automatically add the related files into
% your Matlab's PATH during processing.
%
% Considering different toolboxes may use their own function for the same
% processing which has the same name, the specific method will only be
% loaded just before being used as Matlab does not support function
% overloading. Nebertheless the GUI will automatically do this for you.
%
% Most of the processing methods are provided from other toolboxes,
% inclduing MEDI, STI Suite and FANSI. Please cite the correspoding
% paper(s) when you use this GUI.
%
% **********************************
% * Hierarchy of the GUI structure *
% **********************************
% |---fig
%   |---Tabs
%     |---StepsPanel
%       |---Method(also Edit,Text,etc.)
%         |---Edit, text, Popup, etc.
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 14 September 2017
% Date last modified: 1 June 2018
%
%
function qsm_hub
clear 

% qsm_hub_AddPath
qsm_hub_AddMethodPath;

global h

% set GUI window size
screenSize = get(0,'ScreenSize');
posLeft = round(screenSize(3)/4);
posBottom = round(screenSize(4)/6);
guiSizeHori = round(screenSize(3)/3);
guiSizeVert = round(screenSize(4)*2/3);
if guiSizeHori < 950
    guiSizeHori = 950;
end
if guiSizeVert < 700
    guiSizeVert = 700;
end

% create GUI 
h.fig=figure('Units','pixels','position',[posLeft posBottom guiSizeHori guiSizeVert],...
    'MenuBar','None','Toolbar','None','Name','QSM hub (beta)','NumberTitle','off');

% create Tabs for GUI
h.TabGroup          = uitabgroup(h.fig,'position',[.01 .01 0.98 0.98]);
h.Tabs.QSMHub       = uitab(h.TabGroup,'Title','One-stop QSM processing');
h.Tabs.phaseUnwrap  = uitab(h.TabGroup,'Title','Phase unwrapping');
h.Tabs.bkgRemoval   = uitab(h.TabGroup,'Title','Background field removal');
h.Tabs.qsm          = uitab(h.TabGroup,'Title','QSM');
h.Tabs.utility      = uitab(h.TabGroup,'Title','Utility');

% construct all tabs
%% Phase unwrapping tab
% I/O
h = qsmhub_handle_panel_dataIO(h.Tabs.phaseUnwrap,              h,[0.01 0.8]);
% phase unwrap
h = qsmhub_handle_panel_phaseUnwrap(h.Tabs.phaseUnwrap,         h,[0.01 0.59]);

%% background field removal tab
% I/O
h = qsmhub_handle_panel_dataIO(h.Tabs.bkgRemoval,               h,[0.01 0.8]);
% background field
h = qsmhub_handle_panel_bkgRemoval(h.Tabs.bkgRemoval,           h,[0.01 0.54]);

%% qsm tab
% I/O
h = qsmhub_handle_panel_dataIO(h.Tabs.qsm,                      h,[0.01 0.8]);
% QSM
h = qsmhub_handle_panel_qsm(h.Tabs.qsm,                         h,[0.01 0.54]);

%% utility tab
h = qsmhub_handle_panel_Utility(h.Tabs.utility,                 h,[0.01 0.39]);

%% GUI with QSM one-stop station tab
% I/O
h = qsmhub_handle_panel_dataIO(h.Tabs.QSMHub,                   h,[0.01 0.8]);
% phase unwrap
h = qsmhub_handle_panel_phaseUnwrap(h.Tabs.QSMHub,              h,[0.01 0.59]);
% background field
h = qsmhub_handle_panel_bkgRemoval(h.Tabs.QSMHub,               h,[0.01 0.33]);
% QSM
h = qsmhub_handle_panel_qsm(h.Tabs.QSMHub,                      h,[0.01 0.07]);

% Start button
h.pushbutton_start = uicontrol('Parent',h.Tabs.QSMHub,...
    'Style','pushbutton',...
    'String','Start',...
    'units','normalized','Position',[0.85 0.01 0.1 0.05],...
    'backgroundcolor',get(h.fig,'color'));
% GPU checkbox
h.checkbox_gpu = uicontrol('Parent',h.Tabs.QSMHub,...
    'Style','checkbox',...
    'String','Enable GPU computation',...
    'units','normalized','Position',[0.01 0.01 0.4 0.05],...
    'backgroundcolor',get(h.fig,'color'), 'Enable','off',...
    'TooltipString',['Enable to use GPU for some of the algorithms in qsm_hub. ' ...
                     'Your GPU has to be detectable in Matlab in order to use this feature.']);
if gpuDeviceCount > 0
    set(h.checkbox_gpu, 'Enable', 'on');
end

%% Set Callback functions
set(h.TabGroup,      	'SelectionChangedFcn', {@SwitchTab_Callback})
set(h.pushbutton_start,	'Callback',            {@PushbuttonStart_Callback});

end

%% Callback functions
function SwitchTab_Callback(source,eventdata)
% switch parent handle of StepsPanel based on current tab

global h

switch eventdata.NewValue.Title
    
    % QSM one-stop station tab
    case 'One-stop QSM processing'
        % I/O
        set(h.StepsPanel.dataIO,        'Parent',h.Tabs.QSMHub);
            % This tab supports both DICOM and NIfTI files
            set(h.dataIO.text.input, 'Tooltip',...
                'Input directory contains DICOM (both magnitude and phase files under the same directory) or NIfTI (*phase*.nii* and *magn*.nii*) files');
            % BET is supported with this tab
            set(h.dataIO.checkbox.brainExtraction,'Enable','on');
            % phase invert is supported with this tab
            set(h.dataIO.checkbox.invertPhase,'Enable','on');
        % phase unwrap
        set(h.StepsPanel.phaseUnwrap,   'Parent',h.Tabs.QSMHub,'Position',[0.01 0.59 0.95 0.2]);
        % background field
        set(h.StepsPanel.bkgRemoval,    'Parent',h.Tabs.QSMHub,'Position',[0.01 0.33 0.95 0.25]);
        % QSM
        set(h.StepsPanel.qsm,           'Parent',h.Tabs.QSMHub,'Position',[0.01 0.07 0.95 0.25]);
        % Start pushbutton
        set(h.pushbutton_start,         'Parent',h.Tabs.QSMHub);
        % GPU checkbox
        set(h.checkbox_gpu,             'Parent',h.Tabs.QSMHub);
%         set(h.checkbox_gpu,             'Enable', 'on');
        
    % Phase unwrapping tab
    case 'Phase unwrapping'
        % I/O
        set(h.StepsPanel.dataIO,        'Parent',h.Tabs.phaseUnwrap);
            % This tab supports both DICOM and NIfTI files
            set(h.dataIO.text.input, 'Tooltip',...
                'Input directory contains DICOM (both magnitude and phase files under the same directory) or NIfTI (*phase*.nii* and *magn*.nii*) files');
            % BET is supported with this tab
            set(h.dataIO.checkbox.brainExtraction,'Enable','on');
            % phase invert is supported with this tab
            set(h.dataIO.checkbox.invertPhase,'Enable','on');
        % phase unwrap
        set(h.StepsPanel.phaseUnwrap,   'Parent',h.Tabs.phaseUnwrap,'Position',[0.01 0.59 0.95 0.2]);
        % Start pushbutton
        set(h.pushbutton_start,         'Parent',h.Tabs.phaseUnwrap);
        % GPU checkbox
        set(h.checkbox_gpu,             'Parent',h.Tabs.phaseUnwrap);
%         set(h.checkbox_gpu,             'Enable', 'off');
%         set(h.checkbox_gpu,             'Value',0);
        
    % background field removal tab    
    case 'Background field removal'
        % I/O
        set(h.StepsPanel.dataIO,        'Parent',h.Tabs.bkgRemoval);
            % This tab supports only NIfTI files
            set(h.dataIO.text.input, 'Tooltip',...
                'Input directory contains the unwrapped total field map NIfTI (*totalField*.nii*) files');
            % no BET support with this tab
            set(h.dataIO.checkbox.brainExtraction,  'Enable','off','Value',0);
            set(h.dataIO.edit.maskdir,              'Enable','on');
            set(h.dataIO.button.maskdir,            'Enable','on');
            % phase invert is not supported with this tab
            set(h.dataIO.checkbox.invertPhase,      'Enable','off','Value',0);
        % background field
        set(h.StepsPanel.bkgRemoval,    'Parent',h.Tabs.bkgRemoval,'Position',[0.01 0.54 0.95 0.25]);
        % Start pushbutton
        set(h.pushbutton_start,         'Parent',h.Tabs.bkgRemoval);
        % GPU checkbox
        set(h.checkbox_gpu,             'Parent',h.Tabs.bkgRemoval);
%         set(h.checkbox_gpu,             'Enable', 'on');

    % qsm tab    
    case 'QSM'
        % I/O
        set(h.StepsPanel.dataIO,        'Parent',h.Tabs.qsm);
            % This tab supports only NIfTI files
            set(h.dataIO.text.input, 'Tooltip',...
                'Input directory contains the local field map NIfTI (*localField*.nii*) files; for some QSM methods additional file(s) may also be needed (e.g. *magn*.nii* and *fieldmapSD*.nii*)');
            % no BET support with this tab
            set(h.dataIO.checkbox.brainExtraction,  'Enable','off','Value',0);
            set(h.dataIO.edit.maskdir,              'Enable','on');
            set(h.dataIO.button.maskdir,            'Enable','on');
            % phase invert is not supported with this tab
            set(h.dataIO.checkbox.invertPhase,      'Enable','off','Value',0);
        % QSM
        set(h.StepsPanel.qsm,           'Parent',h.Tabs.qsm,'Position',[0.01 0.54 0.95 0.25]);
        % Start pushbutton
        set(h.pushbutton_start,         'Parent',h.Tabs.qsm);
        % GPU checkbox
        set(h.checkbox_gpu,             'Parent',h.Tabs.qsm);
%         set(h.checkbox_gpu,             'Enable', 'on');
        
end

end

function PushbuttonStart_Callback(source,eventdata)
% determine which step(s) is going to process

global h

% Disable the pushbutton to prevent doubel click
set(source,'Enable','off');

% initialise all possible parameters
subsampling=1;
BFR_tol=1e-4;BFR_depth=4;BFR_peel=2;BFR_iteration=50;
BFR_padSize = 40;        
BFR_radius=4;BFR_alpha=0.01;BFR_threshold=0.03;
QSM_threshold=0.15;QSM_lambda=0.13;QSM_optimise=false;
QSM_tol=1e-3;QSM_maxiter=50;QSM_tol1=0.01;QSM_tol2=0.001;QSM_padsize=[4,4,4];
QSM_mu1=5e-5;QSM_solver='linear';QSM_constraint='tv';QSM_mu2=1;
QSM_radius=5;QSM_zeropad=0;QSM_wData=1;QSM_wGradient=1;QSM_lambdaCSF=100;
QSM_isSMV=false;QSM_merit=false;QSM_isLambdaCSF=false;

% get GPU enable
isGPU = get(h.checkbox_gpu,'Value');

% get I/O GUI input
inputDir        = get(h.dataIO.edit.input,                  'String');
outputBasename  = get(h.dataIO.edit.output,                 'String');
maskFullName    = get(h.dataIO.edit.maskdir,                'String');
isBET           = get(h.dataIO.checkbox.brainExtraction,    'Value');
isInvert        = get(h.dataIO.checkbox.invertPhase,        'Value');

% get phase unwrap GUI input
phaseCombMethod = h.phaseUnwrap.popup.phaseCombMethod.String{h.phaseUnwrap.popup.phaseCombMethod.Value,1};
phaseUnwrap     = h.phaseUnwrap.popup.phaseUnwrap.String{h.phaseUnwrap.popup.phaseUnwrap.Value,1};
isEddyCorrect   = get(h.phaseUnwrap.checkbox.eddyCorrect,'Value');
if get(h.phaseUnwrap.checkbox.excludeMask,'Value')
    excludeMaskThreshold = str2double(get(h.phaseUnwrap.edit.excludeMask,'String'));
else
    excludeMaskThreshold = Inf;
end

% get background field removal GUI input
BFR             = h.bkgRemoval.popup.bkgRemoval.String{h.bkgRemoval.popup.bkgRemoval.Value,1};
refine          = get(h.bkgRemoval.checkbox.bkgRemoval,'Value');

% get QSM GUI input
QSM_method      = h.qsm.popup.qsm.String{h.qsm.popup.qsm.Value,1};

switch phaseUnwrap
    case 'Region growing'
        phaseUnwrap = 'rg';
        
    case 'Graphcut'
        phaseUnwrap = 'gc';
        
    case 'Laplacian STI suite'
        phaseUnwrap = 'laplacian_stisuite';
        
    case '3D best path'
        phaseUnwrap = 'bestpath3d';
        
end

% match the background field removal GUI input to QSMHub input format
% get specific backgroud field removal algorithm parameters
switch BFR
    case 'LBV'
        BFR='lbv';
        BFR_tol         = str2double(get(h.bkgRemoval.LBV.edit.tol,         'String'));
        BFR_depth       = str2double(get(h.bkgRemoval.LBV.edit.depth,       'String'));
        BFR_peel        = str2double(get(h.bkgRemoval.LBV.edit.peel,        'String'));
        
    case 'PDF'
        BFR='pdf';
        BFR_tol         = str2double(get(h.bkgRemoval.PDF.edit.tol,         'String'));
        BFR_iteration   = str2double(get(h.bkgRemoval.PDF.edit.maxIter,     'String'));
        BFR_padSize     = str2double(get(h.bkgRemoval.PDF.edit.padSize,     'String'));

    case 'RESHARP'
        BFR='resharp';
        BFR_radius      = str2double(get(h.bkgRemoval.RESHARP.edit.radius,  'String'));  
        BFR_alpha       = str2double(get(h.bkgRemoval.RESHARP.edit.lambda,  'String'));  
        
    case 'SHARP'
        BFR='sharp';
        BFR_radius      = str2double(get(h.bkgRemoval.SHARP.edit.radius,    'String'));
        BFR_threshold   = str2double(get(h.bkgRemoval.SHARP.edit.threshold, 'String')); 
        
    case 'VSHARP'
        BFR='vsharp';
        maxRadius       = str2double(get(h.bkgRemoval.VSHARP.edit.maxRadius,'String')); 
        minRadius       = str2double(get(h.bkgRemoval.VSHARP.edit.minRadius,'String')); 
        
        BFR_radius      = maxRadius:-1:minRadius;
        
    case 'iHARPERELLA'
        BFR='iharperella';
        BFR_iteration   = str2double(get(h.bkgRemoval.iHARPERELLA.edit.maxIter,'String')); 
        
    case 'VSHARP STI suite'
        BFR='vsharpsti';
        BFR_radius      = str2double(get(h.bkgRemoval.VSHARPSTI.edit.smvSize,'String')); 
        
end

% get QSM algorithm parameters
switch QSM_method
    case 'TKD'
        QSM_method='tkd';
        QSM_threshold   = str2double(get(h.qsm.TKD.edit.threshold,'String'));  
        
    case 'Closed-form solution'
        QSM_method='closedforml2';
        QSM_lambda      = str2double(get(h.qsm.cfs.edit.lambda,'String')); 
        QSM_optimise    = get(h.qsm.cfs.checkbox.lambda,'Value');  
        
    case 'STI suite iLSQR'
        QSM_method='stisuiteilsqr';
        QSM_threshold	= str2double(get(h.qsm.STIiLSQR.edit.threshold,'String')); 
        QSM_maxiter     = str2double(get(h.qsm.STIiLSQR.edit.maxIter,'String'));	
        QSM_tol1        = str2double(get(h.qsm.STIiLSQR.edit.tol1,'String'));  
        QSM_tol2        = str2double(get(h.qsm.STIiLSQR.edit.tol2,'String'));  
        QSM_padsize     = str2double(get(h.qsm.STIiLSQR.edit.padSize,'String'));  
        
        QSM_padsize = [QSM_padsize,QSM_padsize,QSM_padsize];
        
    case 'iLSQR'
        QSM_method='ilsqr';
        QSM_tol         = str2double(get(h.qsm.iLSQR.edit.tol,'String'));        
        QSM_maxiter     = str2double(get(h.qsm.iLSQR.edit.maxIter,'String'));
        QSM_lambda      = str2double(get(h.qsm.iLSQR.edit.lambda,'String'));  
        QSM_optimise    = get(h.qsm.iLSQR.checkbox.lambda,'Value'); 
        
    case 'FANSI'
        QSM_method='fansi';
        QSM_tol         = str2double(get(h.qsm.FANSI.edit.tol,'String'));    
        QSM_lambda      = str2double(get(h.qsm.FANSI.edit.lambda,'String'));   
        QSM_mu1         = str2double(get(h.qsm.FANSI.edit.mu,'String'));    
        QSM_mu2         = str2double(get(h.qsm.FANSI.edit.mu2,'String'));  
        QSM_maxiter     = str2double(get(h.qsm.FANSI.edit.maxIter,'String'));   
        QSM_solver      = h.qsm.FANSI.popup.solver.String{h.qsm.FANSI.popup.solver.Value,1}; 
        QSM_constraint  = h.qsm.FANSI.popup.constraints.String{h.qsm.FANSI.popup.constraints.Value,1}; 
        
    case 'Star-QSM'
        QSM_method='star';
        QSM_padsize     = str2double(get(h.qsm.Star.edit.padSize,'String')); 
        QSM_padsize     = [QSM_padsize,QSM_padsize,QSM_padsize];
        
    case 'MEDI'
        QSM_method='medi_l1';
        QSM_lambda      = str2double(get(h.qsm.MEDI.edit.lambda,'String'));   
        QSM_wData       = str2double(get(h.qsm.MEDI.edit.weightData,'String'));    
        QSM_wGradient   = str2double(get(h.qsm.MEDI.edit.weightGradient,'String')); 
        QSM_zeropad     = str2double(get(h.qsm.MEDI.edit.zeropad,'String'));    
        QSM_radius      = str2double(get(h.qsm.MEDI.edit.smv_radius,'String'));   
        QSM_isSMV       = get(h.qsm.MEDI.checkbox.smv,'Value');   
        QSM_merit       = get(h.qsm.MEDI.checkbox.merit,'Value');     
        QSM_isLambdaCSF = get(h.qsm.MEDI.checkbox.lambda_csf,'Value');        
        QSM_lambdaCSF   = str2double(get(h.qsm.MEDI.edit.lambda_csf,'String'));   
        
end

% run the selected processing step and delay the error (if any)
try 
    % get the parent of current panel, this determine which script is about to run
    switch h.StepsPanel.dataIO.Parent.Title
        case 'One-stop QSM processing'
            % core of QSM one-stop processing
            QSMHub( inputDir,...
                    outputBasename,...
                    maskFullName,...
                    'invert',isInvert,'FSLBet',isBET,'eddy',isEddyCorrect,'GPU',isGPU,...
                    'phase_combine',phaseCombMethod,'unwrap',phaseUnwrap,...
                    'Subsampling',subsampling,'exclude_threshold',excludeMaskThreshold,...
                    'BFR',BFR,'refine',refine,'BFR_tol',BFR_tol,...
                    'depth',BFR_depth,'peel',BFR_peel,'BFR_iteration',BFR_iteration,'BFR_padsize',BFR_padSize,...
                    'BFR_radius',BFR_radius,'BFR_alpha',BFR_alpha,'BFR_threshold',BFR_threshold,...
                    'QSM',QSM_method,'QSM_threshold',QSM_threshold,'QSM_lambda',QSM_lambda,'QSM_optimise',QSM_optimise,...
                    'QSM_tol',QSM_tol,'QSM_iteration',QSM_maxiter,'QSM_tol1',QSM_tol1,'QSM_tol2',QSM_tol2,...
                    'QSM_padsize',QSM_padsize,'QSM_mu',QSM_mu1,'QSM_mu2',QSM_mu2,QSM_solver,QSM_constraint,...
                    'QSM_zeropad',QSM_zeropad,'QSM_wData',QSM_wData,'QSM_wGradient',QSM_wGradient,'QSM_radius',QSM_radius,...
                    'QSM_isSMV',QSM_isSMV,'QSM_merit',QSM_merit,'QSM_isLambdaCSF',QSM_isLambdaCSF,'QSM_lambdaCSF',QSM_lambdaCSF);

        case 'Phase unwrapping'
            % Core of phase unwrapping only 
            UnwrapPhaseMacroIOWrapper(  inputDir,...
                                        outputBasename,...
                                        maskFullName,...
                                        'invert',isInvert,'FSLBet',isBET,'eddy',isEddyCorrect,...
                                        'phase_combine',phaseCombMethod,'unwrap',phaseUnwrap,...
                                        'Subsampling',subsampling,'exclude_threshold',excludeMaskThreshold);

        case 'Background field removal'
            % core of background field removal only
            BackgroundRemovalMacroIOWrapper(inputDir,...
                                            outputBasename,...
                                            maskFullName,...
                                            'GPU',isGPU,...
                                            'BFR',BFR,'refine',refine,'BFR_tol',BFR_tol,...
                                            'depth',BFR_depth,'peel',BFR_peel,'BFR_iteration',BFR_iteration,'BFR_padsize',BFR_padSize,...
                                            'BFR_radius',BFR_radius,'BFR_alpha',BFR_alpha,'BFR_threshold',BFR_threshold);

        case 'QSM'
            % core of QSM only
            QSMMacroIOWrapper(  inputDir,...
                                outputBasename,...
                                maskFullName,...
                                'GPU',isGPU,...
                                'QSM',QSM_method,'QSM_threshold',QSM_threshold,'QSM_lambda',QSM_lambda,'QSM_optimise',QSM_optimise,...
                                'QSM_tol',QSM_tol,'QSM_iteration',QSM_maxiter,'QSM_tol1',QSM_tol1,'QSM_tol2',QSM_tol2,...
                                'QSM_padsize',QSM_padsize,'QSM_mu',QSM_mu1,'QSM_mu2',QSM_mu2,QSM_solver,QSM_constraint,...
                                'QSM_zeropad',QSM_zeropad,'QSM_wData',QSM_wData,'QSM_wGradient',QSM_wGradient,'QSM_radius',QSM_radius,...
                                'QSM_isSMV',QSM_isSMV,'QSM_merit',QSM_merit,'QSM_isLambdaCSF',QSM_isLambdaCSF,'QSM_lambdaCSF',QSM_lambdaCSF);
    end
catch ME
    % re-enable the start button before displaying the error
    set(source,'Enable','on');
    error(ME.message);
end

% get output directory
output_index = strfind(outputBasename, filesep);
outputDir = outputBasename(1:output_index(end));

% generate a log file
GenerateLogFile(h.StepsPanel.dataIO.Parent.Title);

% re-enable the pushbutton
set(source,'Enable','on');

function GenerateLogFile(tab)

switch tab
    case 'One-stop QSM processing'
        fid = fopen([outputDir filesep 'qsm_hub.log'],'w');
        fprintf(fid,'QSMHub(''%s'',...\n',inputDir);
        fprintf(fid,'''%s'',...\n',outputBasename);
        fprintf(fid,'''%s'',...\n',maskFullName);
        fprintf(fid,'''invert'',%i,''FSLBet'',%i,''eddy'',%i,''GPU'',%i,...\n',isInvert,isBET,isEddyCorrect,isGPU);
        fprintf(fid,'''phase_combine'',''%s'',''unwrap'',''%s'',...\n',phaseCombMethod,phaseUnwrap);
        fprintf(fid,'''Subsampling'',%i,''exclude_threshold'',%g,...\n',subsampling,excludeMaskThreshold);
        fprintf(fid,'''BFR'',''%s'',''refine'',%i,''BFR_tol'',%g,...\n',BFR,refine,BFR_tol);
        fprintf(fid,'''depth'',%i,''peel'',%i,''BFR_iteration'',%i,''BFR_padsize'',%i,...\n',BFR_depth,BFR_peel,BFR_iteration,BFR_padSize);
        fprintf(fid,'''BFR_radius'',[%s],''BFR_alpha'',%g,''BFR_threshold'',%g,...\n',num2str(BFR_radius),BFR_alpha,BFR_threshold);
        fprintf(fid,'''QSM'',''%s'',''QSM_threshold'',%g,''QSM_lambda'',%g,''QSM_optimise'',%i,...\n',QSM_method,QSM_threshold,QSM_lambda,QSM_optimise);
        fprintf(fid,'''QSM_tol'',%g,''QSM_iteration'',%i,''QSM_tol1'',%g,''QSM_tol2'',%g,...\n',QSM_tol,QSM_maxiter,QSM_tol1,QSM_tol2);
        fprintf(fid,'''QSM_padsize'',[%s],''QSM_mu'',%g,''QSM_mu2'',%g,''%s'',''%s'',...\n',num2str(QSM_padsize),QSM_mu1,QSM_mu2,QSM_solver,QSM_constraint);
        fprintf(fid,'''QSM_zeropad'',%i,''QSM_wData'',%g,''QSM_wGradient'',%g,''QSM_radius'',%i,...\n',QSM_zeropad,QSM_wData,QSM_wGradient,QSM_radius);
        fprintf(fid,'''QSM_isSMV'',%i,''QSM_merit'',%i,''QSM_isLambdaCSF'',%g,''QSM_lambdaCSF'',%g);\n',QSM_isSMV,QSM_merit,QSM_isLambdaCSF,QSM_lambdaCSF);
        fclose(fid);
        
    case 'Phase unwrapping'
        fid = fopen([outputDir filesep 'qsm_hub.log'],'w');
        fprintf(fid,'UnwrapPhaseMacroIOWrapper(''%s'',...\n',inputDir);
        fprintf(fid,'''%s'',...\n',outputBasename);
        fprintf(fid,'''%s'',...\n',maskFullName);
        fprintf(fid,'''invert'',%i,''FSLBet'',%i,''eddy'',%i,...\n',isInvert,isBET,isEddyCorrect);
        fprintf(fid,'''phase_combine'',''%s'',''unwrap'',''%s'',...\n',phaseCombMethod,phaseUnwrap);
        fprintf(fid,'''Subsampling'',%i,''exclude_threshold'',%g);\n',subsampling,excludeMaskThreshold);
        fclose(fid);
    
    case 'Background field removal'
        fid = fopen([outputDir filesep 'qsm_hub.log'],'w');
        fprintf(fid,'BackgroundRemovalMacroIOWrapper(''%s'',...\n',inputDir);
        fprintf(fid,'''%s'',...\n',outputBasename);
        fprintf(fid,'''%s'',...\n',maskFullName);
        fprintf(fid,'''GPU'',''%i'',...\n',isGPU);
        fprintf(fid,'''BFR'',''%s'',''refine'',%i,''BFR_tol'',%g,...\n',BFR,refine,BFR_tol);
        fprintf(fid,'''depth'',%i,''peel'',%i,''BFR_iteration'',%i,''BFR_padsize'',%i,...\n',BFR_depth,BFR_peel,BFR_iteration,BFR_padSize);
        fprintf(fid,'''BFR_radius'',[%s],''BFR_alpha'',%g,''BFR_threshold'',%g);\n',num2str(BFR_radius),BFR_alpha,BFR_threshold);
        fclose(fid);
    
    case 'QSM'
        fid = fopen([outputDir filesep 'qsm_hub.log'],'w');
        fprintf(fid,'qsmMacroIOWrapper(''%s'',...\n',inputDir);
        fprintf(fid,'''%s'',...\n',outputBasename);
        fprintf(fid,'''%s'',...\n',maskFullName);
        fprintf(fid,'''GPU'',''%i'',...\n',isGPU);
        fprintf(fid,'''QSM'',''%s'',''QSM_threshold'',%g,''QSM_lambda'',%g,''QSM_optimise'',%i,...\n',QSM_method,QSM_threshold,QSM_lambda,QSM_optimise);
        fprintf(fid,'''QSM_tol'',%g,''QSM_iteration'',%i,''QSM_tol1'',%g,''QSM_tol2'',%g,...\n',QSM_tol,QSM_maxiter,QSM_tol1,QSM_tol2);
        fprintf(fid,'''QSM_padsize'',[%s],''QSM_mu'',%g,''QSM_mu2'',%g,''%s'',''%s'',...\n',num2str(QSM_padsize),QSM_mu1,QSM_mu2,QSM_solver,QSM_constraint);
        fprintf(fid,'''QSM_zeropad'',%i,''QSM_wData'',%g,''QSM_wGradient'',%g,''QSM_radius'',%i,...\n',QSM_zeropad,QSM_wData,QSM_wGradient,QSM_radius);
        fprintf(fid,'''QSM_isSMV'',%i,''QSM_merit'',%i,''QSM_isLambdaCSF'',%g,''QSM_lambdaCSF'',%g);\n',QSM_isSMV,QSM_merit,QSM_isLambdaCSF,QSM_lambdaCSF);
        fclose(fid);

end

end

end
