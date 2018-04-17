%% qsm_hub
%
% levels
% ---fig
%  |---Tabs
%    |---StepsPanel
%      |---Edit,Text,etc.

% Description: This is a GUI of QSMHub, which is a pipeline control tool
% for standard QSM processing.
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 10 April 2018
% Date last modified: 17 April 2018
%
%

function qsm_hub
clear 

% qsm_hub_AddPath
qsm_hub_AddMethodPath;

global h

screenSize = get(0,'ScreenSize');
posLeft = round(screenSize(3)/4);
posBottom = round(screenSize(4)/6);
guiSizeHori = round(screenSize(3)/3);
guiSizeVert = round(screenSize(4)*2/3);
if guiSizeHori < 500
    guiSizeHori = 500;
end
if guiSizeVert < 650
    guiSizeVert = 650;
end

fig=figure('Units','pixels','position',[posLeft posBottom guiSizeHori guiSizeVert],...
    'MenuBar','None','Toolbar','None','Name','QSM hub','NumberTitle','off');

h.TabGroup          = uitabgroup(fig,'position',[.01 .01 1 1]);
h.Tabs.QSMHub       = uitab(h.TabGroup,'Title','One-stop QSM processing');
h.Tabs.phaseUnwrap  = uitab(h.TabGroup,'Title','Phase unwrapping');
h.Tabs.bkgRemoval   = uitab(h.TabGroup,'Title','Background field removal');
h.Tabs.qsm          = uitab(h.TabGroup,'Title','QSM');

% initialise all tabs
%% Phase unwrapping tab
% I/O
h = qsmhub_handle_panel_dataIO(h.Tabs.phaseUnwrap,fig,h,[0.01 0.8]);
% phase unwrap
h = qsmhub_handle_panel_phaseUnwrap(h.Tabs.phaseUnwrap,fig,h,[0.01 0.59]);

%% background field removal tab
% I/O
h = qsmhub_handle_panel_dataIO(h.Tabs.bkgRemoval,fig,h,[0.01 0.8]);
% background field
h = qsmhub_handle_panel_bkgRemoval(h.Tabs.bkgRemoval,fig,h,[0.01 0.59]);

%% qsm tab
% I/O
h = qsmhub_handle_panel_dataIO(h.Tabs.qsm,fig,h,[0.01 0.8]);
% QSM
h = qsmhub_handle_panel_qsm(h.Tabs.qsm,fig,h,[0.01 0.59]);

%% GUI with QSM one-stop station tab
% I/O
h = qsmhub_handle_panel_dataIO(h.Tabs.QSMHub,fig,h,[0.01 0.8]);
% phase unwrap
h = qsmhub_handle_panel_phaseUnwrap(h.Tabs.QSMHub,fig,h,[0.01 0.59]);
% background field
h = qsmhub_handle_panel_bkgRemoval(h.Tabs.QSMHub,fig,h,[0.01 0.33]);
% QSM
h = qsmhub_handle_panel_qsm(h.Tabs.QSMHub,fig,h,[0.01 0.07]);
% Start button
h.pushbutton_start = uicontrol('Parent',h.Tabs.QSMHub,'Style','pushbutton',...
    'String','Start',...
    'units','normalized','Position',[0.85 0.01 0.1 0.05],...
    'backgroundcolor',get(fig,'color'));

%% Set Callback
h = SetAllCallbacks(h);
end

%% utils functions
function h=SetAllCallbacks(h)
set(h.TabGroup,                             'SelectionChangedFcn',  {@SwitchTab_Callback})
set(h.dataIO.button.input,                  'Callback',             {@ButtonGetInputDir_Callback});
set(h.dataIO.button.output,                 'Callback',             {@ButtonGetOutputDir_Callback});
set(h.dataIO.checkbox.brainExtraction,      'Callback',             {@CheckboxBrainExtraction_Callback});
set(h.dataIO.button.maskdir,                'Callback',             {@ButtonGetMaskDir_Callback});
set(h.phaseUnwrap.checkbox.excludeMask,     'Callback',             {@CheckboxEditPair_Callback,h.phaseUnwrap.edit.excludeMask,1});
set(h.phaseUnwrap.edit.excludeMask,         'Callback',             {@EditRange01_Callback});
set(h.bkgRemoval.popup.bkgRemoval,      	'Callback',             {@PopupBkgRemoval_Callback});
set(h.bkgRemoval.LBV.edit.depth,            'Callback',             {@EditMinInput_Callback,-1});
set(h.bkgRemoval.LBV.edit.peel,             'Callback',             {@EditMinInput_Callback,0});
set(h.bkgRemoval.LBV.edit.tol,              'Callback',             {@EditMinInput_Callback,0});
set(h.bkgRemoval.PDF.edit.maxIter,        	'Callback',             {@EditMinInput_Callback,0});
set(h.bkgRemoval.PDF.edit.tol,              'Callback',             {@EditMinInput_Callback,0});
set(h.bkgRemoval.PDF.edit.padSize,      	'Callback',             {@EditMinInput_Callback,0});
set(h.bkgRemoval.RESHARP.edit.lambda,     	'Callback',             {@EditMinInput_Callback,0});
set(h.bkgRemoval.RESHARP.edit.radius,     	'Callback',             {@EditMinInput_Callback,0});
set(h.bkgRemoval.SHARP.edit.radius,     	'Callback',             {@EditMinInput_Callback,0});
set(h.bkgRemoval.SHARP.edit.threshold,    	'Callback',             {@EditMinInput_Callback,0});
set(h.bkgRemoval.VSHARP.edit.minRadius,   	'Callback',             {@EditMinInput_Callback,0});
set(h.bkgRemoval.VSHARP.edit.maxRadius,   	'Callback',             {@EditVSHARPMaxRadius_Callback});
set(h.bkgRemoval.VSHARPSTI.edit.smvSize,    'Callback',             {@EditMinInput_Callback,0});
set(h.bkgRemoval.iHARPERELLA.edit.maxIter,	'Callback',             {@EditMinInput_Callback,0});
set(h.qsm.popup.qsm,                        'Callback',             {@PopupQSM_Callback});
set(h.qsm.cfs.checkbox.lambda,              'Callback',             {@CheckboxEditPair_Callback,h.qsm.cfs.edit.lambda,0});
set(h.qsm.iLSQR.checkbox.lambda,            'Callback',             {@CheckboxEditPair_Callback,h.qsm.iLSQR.edit.lambda,0});
set(h.qsm.TKD.edit.threshold,               'Callback',             {@EditRange01_Callback});
set(h.qsm.cfs.edit.lambda,                  'Callback',             {@EditMinInput_Callback,0});
set(h.qsm.Star.edit.padSize,                'Callback',             {@EditMinInput_Callback,0});
set(h.qsm.iLSQR.edit.lambda,                'Callback',             {@EditMinInput_Callback,0});
set(h.qsm.iLSQR.edit.maxIter,               'Callback',             {@EditMinInput_Callback,0});
set(h.qsm.iLSQR.edit.tol,                   'Callback',             {@EditMinInput_Callback,0});
set(h.qsm.STIiLSQR.edit.maxIter,            'Callback',             {@EditMinInput_Callback,0});
set(h.qsm.STIiLSQR.edit.padSize,            'Callback',             {@EditMinInput_Callback,0});
set(h.qsm.STIiLSQR.edit.threshold,          'Callback',             {@EditMinInput_Callback,0});
set(h.qsm.STIiLSQR.edit.tol1,               'Callback',             {@EditMinInput_Callback,0});
set(h.qsm.STIiLSQR.edit.tol2,               'Callback',             {@EditMinInput_Callback,0});
set(h.qsm.FANSI.edit.lambda,                'Callback',             {@EditMinInput_Callback,0});
set(h.qsm.FANSI.edit.maxIter,               'Callback',             {@EditMinInput_Callback,0});
set(h.qsm.FANSI.edit.mu,                    'Callback',             {@EditMinInput_Callback,0});
set(h.qsm.FANSI.edit.tol,                   'Callback',             {@EditMinInput_Callback,0});
set(h.qsm.MEDI.checkbox.smv,                'Callback',             {@CheckboxEditPair_Callback,h.qsm.MEDI.edit.smv_radius,1});
set(h.qsm.MEDI.checkbox.lambda_csf,         'Callback',             {@CheckboxEditPair_Callback,h.qsm.MEDI.edit.lambda_csf,1});
set(h.pushbutton_start,                     'Callback',             {@PushbuttonStart_Callback});
end

%% Callback functions
% common callback functions
function CheckboxEditPair_Callback(source,eventdata,handleToBeDisable,trueValue)
    % get source value
    if source.Value == trueValue
        % if source is equal to trueValue then enables target handle 
        set(handleToBeDisable,'Enable','on');
    else
        % if source do not equal to trueValue then disables target handle 
        set(handleToBeDisable,'Enable','off');
    end

end

function EditMinInput_Callback(source,eventdata,lb)
    if str2double(source.String)<lb
        source.String = num2str(lb);
    end
    
end
function EditRange01_Callback(source,eventdata)
    if str2double(source.String)<0
        source.String = '0';
    end
    if str2double(source.String)>1
        source.String = '1';
    end

end
% specific callback functions
function SwitchTab_Callback(source,eventdata)
global h
% switching parent handle based on current tab
switch eventdata.NewValue.Title
    % QSM one-stop station
    case 'One-stop QSM processing'
        % I/O
        set(h.StepsPanel.dataIO,        'Parent',h.Tabs.QSMHub);
        set(h.dataIO.checkbox.brainExtraction,'Enable','on');
        % phase unwrap
        set(h.StepsPanel.phaseUnwrap,   'Parent',h.Tabs.QSMHub,'Position',[0.01 0.59 0.95 0.2]);
        % background field
        set(h.StepsPanel.bkgRemoval,    'Parent',h.Tabs.QSMHub,'Position',[0.01 0.33 0.95 0.25]);
        % QSM
        set(h.StepsPanel.qsm,           'Parent',h.Tabs.QSMHub,'Position',[0.01 0.07 0.95 0.25]);
        % Start pushbutton
        set(h.pushbutton_start,         'Parent',h.Tabs.QSMHub);
        
    % Phase unwrapping tabs
    case 'Phase unwrapping'
        % I/O
        set(h.StepsPanel.dataIO,        'Parent',h.Tabs.phaseUnwrap);
        set(h.dataIO.checkbox.brainExtraction,'Enable','on');
        % phase unwrap
        set(h.StepsPanel.phaseUnwrap,   'Parent',h.Tabs.phaseUnwrap,'Position',[0.01 0.59 0.95 0.2]);
        % Start pushbutton
        set(h.pushbutton_start,         'Parent',h.Tabs.phaseUnwrap);
        
    % background field removal tab    
    case 'Background field removal'
        % I/O
        set(h.StepsPanel.dataIO,        'Parent',h.Tabs.bkgRemoval);
            % not bet support with this tab
            set(h.dataIO.checkbox.brainExtraction,  'Enable','off','Value',0);
            set(h.dataIO.edit.maskdir,              'Enable','on');
        % background field
        set(h.StepsPanel.bkgRemoval,    'Parent',h.Tabs.bkgRemoval,'Position',[0.01 0.59 0.95 0.25]);
        % Start pushbutton
        set(h.pushbutton_start,         'Parent',h.Tabs.bkgRemoval);

    % qsm tab    
    case 'QSM'
        % I/O
        set(h.StepsPanel.dataIO,        'Parent',h.Tabs.qsm);
            % not bet support with this tab
            set(h.dataIO.checkbox.brainExtraction,  'Enable','off','Value',0);
            set(h.dataIO.edit.maskdir,              'Enable','on');
        % QSM
        set(h.StepsPanel.qsm,           'Parent',h.Tabs.qsm,'Position',[0.01 0.59 0.95 0.25]);
        % Start pushbutton
        set(h.pushbutton_start,         'Parent',h.Tabs.qsm);
        
end

end

function ButtonGetInputDir_Callback(source,eventdata)

global h

pathDir = uigetdir;

if pathDir ~= 0
    set(h.dataIO.edit.input,    'String',pathDir);
    set(h.dataIO.edit.output,   'String',[pathDir filesep 'output']);
end
end

function ButtonGetOutputDir_Callback(source,eventdata)

global h

pathDir = uigetdir;

if pathDir ~= 0
    set(h.dataIO.edit.output,'String',pathDir);
end
end

function CheckboxBrainExtraction_Callback(source,eventdata)

global h

if ~h.dataIO.checkbox.brainExtraction.Value
    set(h.dataIO.button.maskdir,'Enable','on');
    set(h.dataIO.edit.maskdir,  'Enable','on');
else
    set(h.dataIO.button.maskdir,'Enable','off');
    set(h.dataIO.edit.maskdir,  'Enable','off');
end

end
    
function ButtonGetMaskDir_Callback(source,eventdata)

global h

[maskfileName,pathDir] = uigetfile({'*.nii.gz';'*.nii'},'Select mask file');

if pathDir ~= 0
    set(h.dataIO.edit.maskdir,'String',fullfile(pathDir,maskfileName));
end

end

function PopupBkgRemoval_Callback(source,eventdata)

global h

% get selected background removal method
method = source.String{source.Value,1} ;

% switch off all panels
fields = fieldnames(h.bkgRemoval.panel);
for kf = 1:length(fields)
    set(h.bkgRemoval.panel.(fields{kf}),    'Visible','off');
end

% switch on specific panel
switch method
    case 'LBV'
        set(h.bkgRemoval.panel.LBV,         'Visible','on');
        
    case 'PDF'
        set(h.bkgRemoval.panel.PDF,         'Visible','on');

    case 'RESHARP'
        set(h.bkgRemoval.panel.RESHARP,     'Visible','on');

    case 'SHARP'
        set(h.bkgRemoval.panel.SHARP,       'Visible','on');

    case 'VSHARP'
        set(h.bkgRemoval.panel.VSHARP,      'Visible','on');

    case 'VSHARP STI suite'
        set(h.bkgRemoval.panel.VSHARPSTI,   'Visible','on');

    case 'iHARPERELLA'
        set(h.bkgRemoval.panel.iHARPERELLA, 'Visible','on');
end

end

function PopupQSM_Callback(source,eventdata)

global h

% get selected QSM method
method = source.String{source.Value,1} ;

% switch off all panels
fields = fieldnames(h.qsm.panel);
for kf = 1:length(fields)
    set(h.qsm.panel.(fields{kf}),   'Visible','off');
end

% switch on specific panel
switch method
    case 'TKD'
        set(h.qsm.panel.TKD,        'Visible','on');

    case 'Closed-form solution'
        set(h.qsm.panel.cfs,        'Visible','on');

    case 'STI suite iLSQR'
        set(h.qsm.panel.STIiLSQR,   'Visible','on');

    case 'iLSQR'
        set(h.qsm.panel.iLSQR,      'Visible','on');

    case 'FANSI'
        set(h.qsm.panel.FANSI,      'Visible','on');

    case 'Star'
        set(h.qsm.panel.Star,       'Visible','on');

    case 'MEDI'
        set(h.qsm.panel.MEDI,       'Visible','on');
end

end

function EditVSHARPMaxRadius_Callback(source,eventdata)
global h
if str2double(source.String) <= str2double(h.edit_VSHARP_minRadius.String)
    source.String = num2str(str2double(h.edit_VSHARP_minRadius.String) +1);
end
end

function PushbuttonStart_Callback(source,eventdata)

global h

% Disable the pushbutton to prevent dis-kick
set(source,'Enable','off');

% initialise all possible parameters
subsampling=1;
BFR_tol=1e-4;BFR_depth=4;BFR_peel=2;BFR_iteration=50;
% BFR_CGdefault=true;
BFR_padSize = 40;
BFR_radius=4;BFR_alpha=0.01;BFR_threshold=0.03;
QSM_threshold=0.15;QSM_lambda=0.13;QSM_optimise=false;
QSM_tol=1e-3;QSM_maxiter=50;QSM_tol1=0.01;QSM_tol2=0.001;QSM_padsize=[4,4,4];
QSM_mu1=5e-5;QSM_solver='linear';QSM_constraint='tv';
QSM_radius=5;QSM_zeropad=0;QSM_wData=1;QSM_wGradient=1;QSM_lambdaCSF=100;
QSM_isSMV=false;QSM_merit=false;QSM_isLambdaCSF=false;

% get I/O GUI input
inputDir        = get(h.dataIO.edit.input,'String');
outputDir       = get(h.dataIO.edit.output,'String');
maskDir         = get(h.dataIO.edit.maskdir,'String');
isBET           = get(h.dataIO.checkbox.brainExtraction,'Value');

% get phase unwrap GUI input
phaseUnwrap     = h.phaseUnwrap.popup.phaseUnwrap.String{h.phaseUnwrap.popup.phaseUnwrap.Value,1};
isEddyCorrect   = get(h.phaseUnwrap.checkbox.eddyCorrect,'Value');
if get(h.phaseUnwrap.checkbox.excludeMask,'Value')
    excludeMaskThreshold = str2double(get(h.phaseUnwrap.edit.excludeMask,'String'));
else
    excludeMaskThreshold = 1;
end

% get background field removal GUI input
BFR             = h.bkgRemoval.popup.bkgRemoval.String{h.bkgRemoval.popup.bkgRemoval.Value,1};
refine          = get(h.bkgRemoval.checkbox.bkgRemoval,'Value');

% get QSM GUI input
QSM_method      = h.qsm.popup.qsm.String{h.qsm.popup.qsm.Value,1};

% matching the mask GUI input to QSMHub input format
% look for mask file full name
try 
    maskFullName = maskDir;
catch
    maskFullName = [];
end

% matching the phase unwrapping GUI input to QSMHub input format
switch phaseUnwrap
    case 'Region growing'
        phaseUnwrap = 'rg';
        
    case 'Graphcut'
        phaseUnwrap = 'gc';
        
    case 'Laplacian STI suite'
        phaseUnwrap = 'laplacian_stisuite';
        
end

% matching the background field removal GUI input to QSMHub input format
% get specific backgroud field removal algorithm parameters
switch BFR
    case 'LBV'
        BFR='lbv';
        try BFR_tol         = str2double(get(h.bkgRemoval.LBV.edit.tol,'String'));              catch; BFR_tol=1e-4;        end
        try BFR_depth       = str2double(get(h.bkgRemoval.LBV.edit.depth,'String'));            catch; BFR_depth=4;         end
        try BFR_peel        = str2double(get(h.bkgRemoval.LBV.edit.peel,'String'));             catch; BFR_peel=4;          end
        
    case 'PDF'
        BFR='pdf';
        try BFR_tol         = str2double(get(h.bkgRemoval.PDF.edit.tol,'String'));              catch; BFR_tol=1e-2;        end
        try BFR_iteration   = str2double(get(h.bkgRemoval.PDF.edit.maxIter,'String'));          catch; BFR_iteration=50;    end
        try BFR_padSize     = str2double(get(h.bkgRemoval.PDF.edit.padSize,'String'));          catch; BFR_iteration=40;    end
%         try BFR_CGdefault = h.popup_PDF_cgSolver.String{h.popup_PDF_cgSolver.Value,1}; catch; BFR_CGdefault=true; end

    case 'RESHARP'
        BFR='resharp';
        try BFR_radius      = str2double(get(h.bkgRemoval.RESHARP.edit.radius,'String'));       catch; BFR_radius=4;        end
        try BFR_alpha       = str2double(get(h.bkgRemoval.RESHARP.edit.lambda,'String'));       catch; BFR_alpha=0.01;      end
        
    case 'SHARP'
        BFR='sharp';
        try BFR_radius      = str2double(get(h.bkgRemoval.SHARP.edit.radius,'String'));         catch; BFR_radius=4;        end
        try BFR_threshold   = str2double(get(h.bkgRemoval.SHARP.edit.threshold,'String'));      catch; BFR_threshold=0.03;	end
        
    case 'VSHARP'
        BFR='vsharp';
        try maxRadius       = str2double(get(h.bkgRemoval.VSHARP.edit.maxRadius,'String'));     catch; maxRadius=10;        end
        try minRadius       = str2double(get(h.bkgRemoval.VSHARP.edit.minRadius,'String'));     catch; minRadius=3;         end
        
        BFR_radius = maxRadius:-1:minRadius;
        
    case 'iHARPERELLA'
        BFR='iharperella';
        try BFR_iteration   = str2double(get(h.bkgRemoval.iHARPERELLA.edit.maxIter,'String'));  catch; BFR_iteration=100;   end 
        
    case 'VSHARP STI suite'
        BFR='vsharpsti';
        try BFR_radius      = str2double(get(h.bkgRemoval.VSHARPSTI.edit.smvSize,'String'));    catch; BFR_radius=12;       end
        
end

% get QSM algorithm parameters
switch QSM_method
    case 'TKD'
        QSM_method='tkd';
        try QSM_threshold   = str2double(get(h.qsm.TKD.edit.threshold,'String'));       catch; QSM_threshold=0.15;  end
        
    case 'Closed-form solution'
        QSM_method='closedforml2';
        try QSM_lambda      = str2double(get(h.qsm.cfs.edit.lambda,'String'));          catch; QSM_lambda=0.13;     end
        try QSM_optimise    = get(h.qsm.cfs.edit.lambda,'Value');                       catch; QSM_optimise=false;  end
        
    case 'STI suite iLSQR'
        QSM_method='stisuiteilsqr';
        try QSM_threshold	= str2double(get(h.qsm.STIiLSQR.edit.threshold,'String'));  catch; QSM_threshold=0.01;  end
        try QSM_maxiter     = str2double(get(h.qsm.STIiLSQR.edit.maxIter,'String'));	catch; QSM_maxiter=100;     end
        try QSM_tol1        = str2double(get(h.qsm.STIiLSQR.edit.tol1,'String'));       catch; QSM_tol1=0.01;       end
        try QSM_tol2        = str2double(get(h.qsm.STIiLSQR.edit.tol2,'String'));       catch; QSM_tol2=0.001;      end
        try QSM_padsize     = str2double(get(h.qsm.STIiLSQR.edit.padSize,'String'));    catch; QSM_padsize=4;       end
        
        QSM_padsize = [QSM_padsize,QSM_padsize,QSM_padsize];
        
    case 'iLSQR'
        QSM_method='ilsqr';
        try QSM_tol         = str2double(get(h.qsm.iLSQR.edit.tol,'String'));           catch; QSM_tol=0.001;       end
        try QSM_maxiter     = str2double(get(h.qsm.iLSQR.edit.maxIter,'String'));       catch; QSM_maxiter=100;     end
        try QSM_lambda      = str2double(get(h.qsm.iLSQR.edit.lambda,'String'));        catch; QSM_lambda=0.13;     end
        try QSM_optimise    = get(h.qsm.iLSQR.edit.lambda,'Value');                     catch; QSM_optimise=false;  end 
        
    case 'FANSI'
        QSM_method='fansi';
        try QSM_tol         = str2double(get(h.qsm.FANSI.edit.tol,'String'));           catch; QSM_tol=1;           end
        try QSM_lambda      = str2double(get(h.qsm.FANSI.edit.lambda,'String'));        catch; QSM_lambda=3e-5;     end
        try QSM_mu1         = str2double(get(h.qsm.FANSI.edit.mu,'String'));            catch; QSM_mu1=5e-5;        end
        try QSM_maxiter     = str2double(get(h.qsm.FANSI.edit.maxIter,'String'));       catch; QSM_maxiter=50;      end
        try 
            QSM_solver      = h.qsm.FANSI.popup.solver.String{h.qsm.FANSI.popup.solver.Value,1}; 
        catch
            QSM_solver      ='linear'; 
        end 
        try 
            QSM_constraint  = h.qsm.FANSI.popup.constraints.String{h.qsm.FANSI.popup.constraints.Value,1}; 
        catch
            QSM_constraint  ='tv'; 
        end 
        
    case 'Star'
        QSM_method='star';
        try QSM_threshold   = str2double(get(h.qsm.Star.edit.padSize,'String'));        catch; QSM_padsize=4;       end
        
    case 'MEDI'
        QSM_method='medi_l1';
        try QSM_lambda      = str2double(get(h.qsm.MEDI.edit.lambda,'String'));         catch; QSM_lambda=1000;     end 
        try QSM_wData       = str2double(get(h.qsm.MEDI.edit.weightData,'String'));     catch; QSM_wData=1;         end 
        try QSM_wGradient   = str2double(get(h.qsm.MEDI.edit.weightGradient,'String')); catch; QSM_wGradient=1;     end 
        try QSM_zeropad     = str2double(get(h.qsm.MEDI.edit.zeropad,'String'));        catch; QSM_zeropad=0;       end 
        try QSM_radius      = str2double(get(h.qsm.MEDI.edit.smv_radius,'String'));     catch; QSM_radius=5;        end 
        try QSM_isSMV       = get(h.qsm.MEDI.checkbox.smv,'Value');                     catch; QSM_isSMV=0;         end 
        try QSM_merit       = get(h.qsm.MEDI.checkbox.merit,'Value');                   catch; QSM_merit=0;         end 
        try QSM_isLambdaCSF = get(h.qsm.MEDI.checkbox.lambda_csf,'Value');              catch; QSM_isLambdaCSF=0;   end 
        try QSM_lambdaCSF   = str2double(get(h.qsm.MEDI.edit.lambda_csf,'String'));     catch; QSM_lambdaCSF=100;   end 
        
end

% get the parent of current panel, this determine which script is about to
% run
switch h.StepsPanel.dataIO.Parent.Title
    case 'One-stop QSM processing'
        % core for QSM one stop processing
        QSMHub(inputDir,outputDir,'FSLBet',isBET,'mask',maskFullName,'unwrap',phaseUnwrap,...
            'Subsampling',subsampling,'BFR',BFR,'refine',refine,'BFR_tol',BFR_tol,...
            'depth',BFR_depth,'peel',BFR_peel,'BFR_iteration',BFR_iteration,'BFR_padsize',BFR_padSize,...
            'BFR_radius',BFR_radius,'BFR_alpha',BFR_alpha,'BFR_threshold',BFR_threshold,...
            'QSM',QSM_method,'QSM_threshold',QSM_threshold,'QSM_lambda',QSM_lambda,'QSM_optimise',QSM_optimise,...
            'QSM_tol',QSM_tol,'QSM_iteration',QSM_maxiter,'QSM_tol1',QSM_tol1,'QSM_tol2',QSM_tol2,...
            'QSM_padsize',QSM_padsize,'QSM_mu',QSM_mu1,QSM_solver,QSM_constraint,'exclude_threshold',excludeMaskThreshold,...
            'QSM_zeropad',QSM_zeropad,'QSM_wData',QSM_wData,'QSM_wGradient',QSM_wGradient,'QSM_radius',QSM_radius,...
            'QSM_isSMV',QSM_isSMV,'QSM_merit',QSM_merit,'QSM_isLambdaCSF',QSM_isLambdaCSF,'QSM_lambdaCSF',QSM_lambdaCSF,'eddy',isEddyCorrect);
    
    case 'Phase unwrapping'
        UnwrapPhaseIOMacro(inputDir,outputDir,'FSLBet',isBET,'mask',maskFullName,'unwrap',phaseUnwrap,...
            'Subsampling',subsampling,'exclude_threshold',excludeMaskThreshold,'eddy',isEddyCorrect);
    
    case 'Background field removal'
        BackgroundRemovalIOMacro(inputDir,outputDir,'mask',maskFullName,...
            'BFR',BFR,'refine',refine,'BFR_tol',BFR_tol,...
            'depth',BFR_depth,'peel',BFR_peel,'BFR_iteration',BFR_iteration,'BFR_padsize',BFR_padSize,...
            'BFR_radius',BFR_radius,'BFR_alpha',BFR_alpha,'BFR_threshold',BFR_threshold);
    
    case 'QSM'
        QSMIOMacro(inputDir,outputDir,'mask',maskFullName,...
            'QSM',QSM_method,'QSM_threshold',QSM_threshold,'QSM_lambda',QSM_lambda,'QSM_optimise',QSM_optimise,...
            'QSM_tol',QSM_tol,'QSM_iteration',QSM_maxiter,'QSM_tol1',QSM_tol1,'QSM_tol2',QSM_tol2,...
            'QSM_padsize',QSM_padsize,'QSM_mu',QSM_mu1,QSM_solver,QSM_constraint,...
            'QSM_zeropad',QSM_zeropad,'QSM_wData',QSM_wData,'QSM_wGradient',QSM_wGradient,'QSM_radius',QSM_radius,...
            'QSM_isSMV',QSM_isSMV,'QSM_merit',QSM_merit,'QSM_isLambdaCSF',QSM_isLambdaCSF,'QSM_lambdaCSF',QSM_lambdaCSF);
end

% re-enable the pushbutton
set(source,'Enable','on');

end
