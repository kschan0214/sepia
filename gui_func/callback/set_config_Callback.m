
function set_config_Callback(config_filename,h)

isModifyGeneral = false;
isModifyUnwrap  = false;
isModifyBFR     = false;
isModifyQSM     = false;

A = fileread(config_filename);

str_pattern	= '.general.';
idx         = regexp(A,str_pattern,'once');
if ~isempty(idx)
    isModifyGeneral = true;
end
str_pattern	= '.unwrap.';
idx         = regexp(A,str_pattern,'once');
if ~isempty(idx)
    isModifyUnwrap = true;
end
str_pattern	= '.bfr.';
idx         = regexp(A,str_pattern,'once');
if ~isempty(idx)
    isModifyBFR = true;
end
str_pattern	= '.qsm.';
idx         = regexp(A,str_pattern,'once');
if ~isempty(idx)
    isModifyQSM = true;
end 

%% Panel: I/O
if isModifyGeneral
    
% invert phase data
str_pattern	= '.general.isInvert';
val      	= get_num_as_string(A, str_pattern,'=',';');
set_non_nan_value(h.dataIO.checkbox.invertPhase, 'Value', str2double(val));

% BET
str_pattern	= '.general.isBET';
val      	= get_num_as_string(A, str_pattern,'=',';');
set_non_nan_value(h.dataIO.checkbox.brainExtraction, 'Value', str2double(val));
% trigger checkout callback
feval(h.dataIO.checkbox.brainExtraction.Callback{1},h.dataIO.checkbox.brainExtraction,[],h);
if str2double(val)
    % -f
    str_pattern	= '.general.fractional_threshold';
    val      	= get_num_as_string(A, str_pattern,'=',';');
    set_non_nan_value(h.dataIO.edit.fractionalThres, 'String', val);
    
    % -g
    str_pattern	= '.general.gradient_threshold';
    val      	= get_num_as_string(A, str_pattern,'=',';');
    set_non_nan_value(h.dataIO.edit.gradientThres, 'String', val);
end

end

%% Panel: Total field
if isModifyUnwrap
    
% echo combination
str_pattern	= '.unwrap.echoCombMethod';
val	= get_string_as_string(A, str_pattern);
if ~isnan(val)
    switch lower(val)
        case 'optimum weights'
            val = 1;
        case 'medi nonlinear fit'
            val = 2;
        case 'medi nonlinear fit (bipolar)'
            val = 3;
    end
    
    % change method popup manu
    set_non_nan_value(h.phaseUnwrap.popup.phaseCombMethod,'Value',val);
    % trigger popup callback to switch method panel
    feval(h.phaseUnwrap.popup.phaseCombMethod.Callback{1},h.phaseUnwrap.popup.phaseCombMethod,[],h)
end

% phase unwrap
str_pattern	= '.unwrap.unwrapMethod';
val	= get_string_as_string(A, str_pattern);
if ~isnan(val)
    switch lower(val)
        case 'laplacian'
            val = 1;
        case 'laplacian_stisuite'
            val = 2;
        case 'bestpath3d'
            val = 3;
        case 'rg'
            val = 4;
        case 'gc'
            val = 5;
        case 'segue'
            val = 6;
    end
    
    % change method popup manu
    set_non_nan_value(h.phaseUnwrap.popup.phaseUnwrap,'Value',val);
    % trigger popup callback to switch method panel
    feval(h.phaseUnwrap.popup.phaseUnwrap.Callback{1},h.phaseUnwrap.popup.phaseUnwrap,[],h)
end

% eddy corr
str_pattern	= '.unwrap.isEddyCorrect';
val      	= get_num_as_string(A, str_pattern,'=',';');
set_non_nan_value(h.phaseUnwrap.checkbox.eddyCorrect, 'Value', str2double(val));

% exlusion threshold
str_pattern	= '.unwrap.excludeMaskThreshold';
val      	= get_num_as_string(A, str_pattern,'=',';');
if isinf(str2double(val))
    % trigger checkbox
    set(h.phaseUnwrap.checkbox.excludeMask,'Value',0);
    feval(h.phaseUnwrap.checkbox.excludeMask.Callback{1},h.phaseUnwrap.checkbox.excludeMask,[],{h.phaseUnwrap.edit.excludeMask,h.phaseUnwrap.popup.excludeMethod},1);
else
    % trigger checkbox
    set(h.phaseUnwrap.checkbox.excludeMask,'Value',1);
    feval(h.phaseUnwrap.checkbox.excludeMask.Callback{1},h.phaseUnwrap.checkbox.excludeMask,[],{h.phaseUnwrap.edit.excludeMask,h.phaseUnwrap.popup.excludeMethod},1);
    
    % modifiy edit field value
    set_non_nan_value(h.phaseUnwrap.edit.excludeMask, 'String', val);
    
    % popup manu for thresholding method
    str_pattern	= '.unwrap.excludeMethod';
    val      	= get_string_as_string(A, str_pattern);
    if ~isnan(val)
        switch lower(val)
            case 'weighting map'
                val = 1;
            case 'brain mask'
                val = 2;
        end
        % change method popup manu
        set_non_nan_value(h.phaseUnwrap.popup.excludeMethod,'Value',val);
    end
end

% save unwrap
str_pattern	= '.unwrap.isSaveUnwrappedEcho';
val      	= get_num_as_string(A, str_pattern,'=',';');
set_non_nan_value(h.phaseUnwrap.checkbox.saveEchoPhase, 'Value', str2double(val));

end

%% Panel: BFR
if isModifyBFR
    
% B1 polyfit
str_pattern	= '.bfr.refine';
val      	= get_num_as_string(A, str_pattern,'=',';');
set_non_nan_value(h.bkgRemoval.checkbox.bkgRemoval, 'Value', str2double(val));

% Erosion radius
str_pattern	= '.bfr.erode_radius';
val      	= get_num_as_string(A, str_pattern,'=',';');
set_non_nan_value(h.bkgRemoval.edit.imerode, 'String', val);

% method
str_pattern	= '.bfr.method';
BFR_method	= get_string_as_string(A, str_pattern);
if ~isnan(BFR_method)
    switch lower(BFR_method)
        case 'lbv'
            set_bfr_lbv(A,h);           BFR_method = 1;
        case 'pdf'
            set_bfr_pdf(A,h);           BFR_method = 2;
        case 'resharp'
            set_bfr_resharp(A,h);    	BFR_method = 3;
        case 'sharp'
            set_bfr_sharp(A,h);         BFR_method = 4;
        case 'vsharpsti'
            set_bfr_vsharpsti(A,h);    	BFR_method = 5;
        case 'vsharp'
            set_bfr_vsharp(A,h);    	BFR_method = 6;
        case 'iharperella'
            set_bfr_iharperella(A,h); 	BFR_method = 7;

    end
    
    % change method popup manu
    set_non_nan_value(h.bkgRemoval.popup.bkgRemoval,'Value',BFR_method);
    % trigger popup callback to switch method panel
    feval(h.bkgRemoval.popup.bkgRemoval.Callback{1},h.bkgRemoval.popup.bkgRemoval,[],h)
end

end
%% Panel: QSM
if isModifyQSM
    
% reference tissue
str_pattern	= '.qsm.reference_tissue';
val      	= get_string_as_string(A, str_pattern);
switch lower(val)
    case 'none'
        set_non_nan_value(h.qsm.popup.tissue,'Value',1)
    case 'brain mask'
        set_non_nan_value(h.qsm.popup.tissue,'Value',2)
    case 'csf'
        set_non_nan_value(h.qsm.popup.tissue,'Value',3)
end
% method
str_pattern	= '.qsm.method';
QSM_method	= get_string_as_string(A, str_pattern);
if ~isnan(QSM_method)
    switch lower(QSM_method)
        case 'tkd'
            set_qsm_tkd(A,h);           QSM_method = 1;
        case 'closedforml2'
            set_qsm_cfs(A,h);           QSM_method = 2;
        case 'ndi'
            set_qsm_ndi(A,h);           QSM_method = 3;
        case 'stisuiteilsqr'
            set_qsm_stisuiteilsqr(A,h); QSM_method = 4;
        case 'ilsqr'
            set_qsm_ilsqr(A,h);         QSM_method = 5;
        case 'fansi'
            set_qsm_fansi(A,h);         QSM_method = 6;
        case 'star'
            set_qsm_starqsm(A,h);       QSM_method = 7;
        case 'medi_l1'
            set_qsm_medi(A,h);          QSM_method = 8;

    end
    
    % change method popup manu
    set_non_nan_value(h.qsm.popup.qsm,'Value',QSM_method);
    % trigger popup callback to switch method panel
    feval(h.qsm.popup.qsm.Callback{1},h.qsm.popup.qsm,[],h)
end

end

end

%% get value from config file
function str = get_num_as_string(A, str_pattern, start_indicator, end_indicator)

% get the last position of thw string
str_end_idx = regexp(A,str_pattern,'end');

% check if there is more than one parameter with similar pattern
if length(str_end_idx) > 1
    for k = 1:length(str_end_idx)
        next_character = A(str_end_idx(k)+1);
        if isspace(next_character) || strcmp(next_character,'=') % if it is space to equal symbol
            str_end_idx = str_end_idx(k);
            break
        end
    end
end

if ~isempty(str_end_idx)

    indicatorS_idx  = regexp(A,start_indicator);
    indicatorE_idx 	= regexp(A,end_indicator);

    % get all characters between the indicators
    str = A(indicatorS_idx(find(indicatorS_idx > str_end_idx, 1 ))+1:indicatorE_idx(find(indicatorE_idx > str_end_idx, 1 ))-1);

    % remove all white space
    str = str(find(~isspace(str)));
else
    str = nan;
end

end

function str = get_string_as_string(A, str_pattern)

% get the last position of thw string
str_end_idx = regexp(A,str_pattern,'end');

if ~isempty(str_end_idx)

    indicator_idx  = regexp(A,'''');

    % get all characters between the indicators
    str = A(indicator_idx(find(indicator_idx > str_end_idx, 1 ))+1:indicator_idx(find(indicator_idx > str_end_idx, 1 )+1)-1);
else
    str = nan;
end

end

% only set the parameter when it is valid
function set_non_nan_value(h,field,val)
    if ~isnan(val)
        set(h,field,val);
    end
end

%% LBV
function set_bfr_lbv(A,h)

str_pattern = {'.bfr.tol',...
               '.bfr.depth',...
               '.bfr.peel'};

action_handle = {h.bkgRemoval.LBV.edit.tol,...
                 h.bkgRemoval.LBV.edit.depth,...
                 h.bkgRemoval.LBV.edit.peel};
           
% first 3 edit fields
for k = 1:length(action_handle)
    pattern_curr    = str_pattern{k};
    val             = get_num_as_string(A, pattern_curr, '=', ';');
    set_non_nan_value(action_handle{k},'String',val)
end

end

%% PDF
function set_bfr_pdf(A,h)

str_pattern = {'.bfr.tol',...
               '.bfr.iteration',...
               '.bfr.padSize'};

action_handle = {h.bkgRemoval.PDF.edit.tol,...
                 h.bkgRemoval.PDF.edit.maxIter,...
                 h.bkgRemoval.PDF.edit.padSize };
           
% first 3 edit fields
for k = 1:length(action_handle)
    pattern_curr    = str_pattern{k};
    val             = get_num_as_string(A, pattern_curr, '=', ';');
    set_non_nan_value(action_handle{k},'String',val)
end

end

%% RESHARP
function set_bfr_resharp(A,h)

str_pattern = {'.bfr.radius',...
               '.bfr.alpha'};

action_handle = {h.bkgRemoval.RESHARP.edit.radius,...
                 h.bkgRemoval.RESHARP.edit.lambda};
           
for k = 1:length(action_handle)
    pattern_curr    = str_pattern{k};
    val             = get_num_as_string(A, pattern_curr, '=', ';');
    set_non_nan_value(action_handle{k},'String',val)
end

end

%% SHARP
function set_bfr_sharp(A,h)

str_pattern = {'.bfr.radius',...
               '.bfr.threshold'};

action_handle = {h.bkgRemoval.SHARP.edit.radius,...
                 h.bkgRemoval.SHARP.edit.threshold};
           
for k = 1:length(action_handle)
    pattern_curr    = str_pattern{k};
    val             = get_num_as_string(A, pattern_curr, '=', ';');
    set_non_nan_value(action_handle{k},'String',val)
end

end

%% VSHARP STI suite
function set_bfr_vsharpsti(A,h)

str_pattern = {'.bfr.radius'};

action_handle = {h.bkgRemoval.VSHARPSTI.edit.smvSize};
           
for k = 1:length(action_handle)
    pattern_curr    = str_pattern{k};
    val             = get_num_as_string(A, pattern_curr, '=', ';');
    set_non_nan_value(action_handle{k},'String',val)
end

end

%% VSHARP 
function set_bfr_vsharp(A,h)

str_pattern = {'.bfr.radius'};

action_handle = {h.bkgRemoval.VSHARP.edit.maxRadius,...
                 h.bkgRemoval.VSHARP.edit.minRadius};
           
% max radius
pattern_curr    = str_pattern{1};
val             = get_num_as_string(A, pattern_curr, '[', ':');
set_non_nan_value(action_handle{1},'String',val)

% min radius
val             = get_num_as_string(A, pattern_curr, ':', ']');
idx             = regexp(val,':');
val             = val(idx+1:end);
set_non_nan_value(action_handle{2},'String',val)


end

%% iharperella
function set_bfr_iharperella(A,h)

str_pattern = {'.bfr.iteration'};

action_handle = {h.bkgRemoval.iHARPERELLA.edit.maxIter};
           
for k = 1:length(action_handle)
    pattern_curr    = str_pattern{k};
    val             = get_num_as_string(A, pattern_curr, '=', ';');
    set_non_nan_value(action_handle{k},'String',val)
end

end

%% TKD
function set_qsm_tkd(A,h)

str_pattern = {'.qsm.threshold'};

action_handle = {h.qsm.TKD.edit.threshold};
           
k = 1;
pattern_curr    = str_pattern{k};
val             = get_num_as_string(A, pattern_curr, '=', ';');
set_non_nan_value(action_handle{k},'String',val);

end

%% CFS
function set_qsm_cfs(A,h)

str_pattern = {'.qsm.lambda',...
               '.qsm.optimise'};

action_handle = {h.qsm.cfs.edit.lambda,...
                 h.qsm.cfs.checkbox.lambda};
           
% Lambda
k = 1;
pattern_curr    = str_pattern{k};
val             = get_num_as_string(A, pattern_curr, '=', ';');
set_non_nan_value(action_handle{k},'String',val)

% L-curve optimisation
k = k+1;
pattern_curr    = str_pattern{k};
val             = get_num_as_string(A, pattern_curr, '=', ';');
set_non_nan_value(action_handle{k}, 'Value', str2double(val))
% trigger popup callback to switch method panel
feval(action_handle{k}.Callback{1},action_handle{k},[],action_handle{k-1},0);

end

%% NDI
function set_qsm_ndi(A,h)

str_pattern = {'.qsm.tol',...
               '.qsm.maxiter',...
               '.qsm.stepSize'};

action_handle = {h.qsm.NDI.edit.tol,...
                 h.qsm.NDI.edit.maxIter,...
                 h.qsm.NDI.edit.stepSize};
           
% first 3 edit fields
for k = 1:length(action_handle)
    pattern_curr    = str_pattern{k};
    val             = get_num_as_string(A, pattern_curr, '=', ';');
    set_non_nan_value(action_handle{k},'String',val)
end

end

%% iLSQR (STI suite)
function set_qsm_stisuiteilsqr(A,h)

str_pattern = {'.qsm.threshold',...
               '.qsm.maxiter',...
               '.qsm.tol1',...
               '.qsm.tol2',...
               '.qsm.padsize'};

action_handle = {h.qsm.STIiLSQR.edit.threshold,...
                 h.qsm.STIiLSQR.edit.maxIter,...
                 h.qsm.STIiLSQR.edit.tol1,...
                 h.qsm.STIiLSQR.edit.tol2,...
                 h.qsm.STIiLSQR.edit.padSize};
           
for k = 1:length(action_handle)-1
    pattern_curr    = str_pattern{k};
    val             = get_num_as_string(A, pattern_curr, '=', ';');
    set_non_nan_value(action_handle{k},'String',val)
end

k = k+1;
pattern_curr    = str_pattern{k};
val             = get_num_as_string(A, pattern_curr, '*', ';');
set_non_nan_value(action_handle{k},'String',val)

end

%% iLSQR
function set_qsm_ilsqr(A,h)

str_pattern = {'.qsm.tol',...
               '.qsm.maxiter',...
               '.qsm.lambda',...
               '.qsm.optimise'};

action_handle = {h.qsm.iLSQR.edit.tol,...
                 h.qsm.iLSQR.edit.maxIter,...
                 h.qsm.iLSQR.edit.lambda,...
                 h.qsm.iLSQR.checkbox.lambda};
           
% first 3 edit fields
for k = 1:3
    pattern_curr    = str_pattern{k};
    val             = get_num_as_string(A, pattern_curr, '=', ';');
    set_non_nan_value(action_handle{k},'String',val)
end

% L-curve optimisation
k = k+1;
pattern_curr    = str_pattern{k};
val             = get_num_as_string(A, pattern_curr, '=', ';');
set_non_nan_value(action_handle{k}, 'Value', str2double(val))
% trigger popup callback to switch method panel
feval(action_handle{k}.Callback{1},action_handle{k},[],action_handle{k-1},0);

% if action_handle{k}.Value
%     set(action_handle{k-1},'Enable','on')
% else
%     set(action_handle{k-1},'Enable','off');
% end

end

%% FANSI
function set_qsm_fansi(A,h)

str_pattern = {'.qsm.tol',...
               '.qsm.maxiter',...
               '.qsm.lambda',...
               '.qsm.mu1',...
               '.qsm.mu2',...
               '.qsm.solver',...
               '.qsm.constraint',...
               '.qsm.gradient_mode',...
               '.qsm.isWeakHarmonic',...
               '.qsm.beta',...
               '.qsm.muh'};

action_handle = {h.qsm.FANSI.edit.tol,...
                 h.qsm.FANSI.edit.maxIter,...
                 h.qsm.FANSI.edit.lambda,...
                 h.qsm.FANSI.edit.mu,...
                 h.qsm.FANSI.edit.mu2,...
                 h.qsm.FANSI.popup.solver,...
                 h.qsm.FANSI.popup.constraints,...
                 h.qsm.FANSI.popup.gradientMode,...
                 h.qsm.FANSI.checkbox.isWeakHarmonic,...
                 h.qsm.FANSI.edit.beta,...
                 h.qsm.FANSI.edit.muh};
           
% first 5 edit fields
for k = 1:5
    pattern_curr    = str_pattern{k};
    val             = get_num_as_string(A, pattern_curr, '=', ';');
    set_non_nan_value(action_handle{k},'String',val)
end

% solver
k = k+1;
pattern_curr    = str_pattern{k};
val             = get_string_as_string(A, pattern_curr);
switch lower(val)
    case 'linear'
        set_non_nan_value(action_handle{k},'Value',1)
    case 'non-linear'
        set_non_nan_value(action_handle{k},'Value',2)
end

% constraint
k = k+1;
pattern_curr    = str_pattern{k};
val             = get_string_as_string(A, pattern_curr);
switch lower(val)
    case 'tv'
        set_non_nan_value(action_handle{k},'Value',1)
    case 'tgv'
        set_non_nan_value(action_handle{k},'Value',2)
end

% gradient mode
k = k+1;
pattern_curr    = str_pattern{k};
val             = get_string_as_string(A, pattern_curr);
switch lower(val)
    case 'vector field'
        set_non_nan_value(action_handle{k},'Value',1)
    case 'l1 norm'
        set_non_nan_value(action_handle{k},'Value',2)
    case 'l2 norm'
        set_non_nan_value(action_handle{k},'Value',3)
end

% weak field harmonic
k = k+1;
pattern_curr    = str_pattern{k};
val             = get_num_as_string(A, pattern_curr, '=', ';');
set_non_nan_value(action_handle{k}, 'Value', str2double(val))
% trigger popup callback to switch method panel
feval(action_handle{k}.Callback{1},action_handle{k},[],{action_handle{k+1},action_handle{k+2}},1);

for k = 10:11
    pattern_curr    = str_pattern{k};
    val             = get_num_as_string(A, pattern_curr, '=', ';');
    set_non_nan_value(action_handle{k},'String',val);
end

end

%% Star-QSM
function set_qsm_starqsm(A,h)

str_pattern = {'.qsm.padsize'};

action_handle = {h.qsm.Star.edit.padSize};
          
k = 1;
pattern_curr    = str_pattern{k};
val             = get_num_as_string(A, pattern_curr, '*', ';');
set_non_nan_value(action_handle{k},'String',val)

end

%% MEDI
function set_qsm_medi(A,h)

str_pattern = {'.qsm.lambda',...
               '.qsm.wData',...
               '.qsm.zeropad',...
               '.qsm.percentage',...
               '.qsm.isSMV',...
               '.qsm.radius',...
               '.qsm.merit',...
               '.qsm.isLambdaCSF',...
               '.qsm.lambdaCSF'};

action_handle = {h.qsm.MEDI.edit.lambda,...
                 h.qsm.MEDI.edit.weightData,...
                 h.qsm.MEDI.edit.zeropad,...
                 h.qsm.MEDI.edit.percentage,...
                 h.qsm.MEDI.checkbox.smv,...
                 h.qsm.MEDI.edit.smv_radius,...
                 h.qsm.MEDI.checkbox.merit,...
                 h.qsm.MEDI.checkbox.lambda_csf,...
                 h.qsm.MEDI.edit.lambda_csf};
           

for k = 1:4
    pattern_curr    = str_pattern{k};
    val             = get_num_as_string(A, pattern_curr, '=', ';');
    set_non_nan_value(action_handle{k},'String',val)
end

% SMV
k = k+1;
pattern_curr    = str_pattern{k};
val             = get_num_as_string(A, pattern_curr, '=', ';');
set_non_nan_value(action_handle{k}, 'Value', str2double(val))
% trigger popup callback to switch method panel
feval(action_handle{k}.Callback{1},action_handle{k},[],action_handle{k+1},1);

% SMV radius
k = k+1;
pattern_curr    = str_pattern{k};
val             = get_num_as_string(A, pattern_curr, '=', ';');
set_non_nan_value(action_handle{k},'String',val);


% Merit
k = k+1;
pattern_curr    = str_pattern{k};
val             = get_num_as_string(A, pattern_curr, '=', ';');
set_non_nan_value(action_handle{k}, 'Value', str2double(val))

% lambda CSF
k = k+1;
pattern_curr    = str_pattern{k};
val             = get_num_as_string(A, pattern_curr, '=', ';');
set_non_nan_value(action_handle{k}, 'Value', str2double(val))
feval(action_handle{k}.Callback{1},action_handle{k},[],action_handle{k+1},1);

% lambda CSF value
k = k+1;
pattern_curr    = str_pattern{k};
val             = get_num_as_string(A, pattern_curr, '=', ';');
set_non_nan_value(action_handle{k},'String',val);


end