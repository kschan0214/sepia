%% function output = function_name(input)
%
% Usage:
%
% Input
% --------------
%
% Output
% --------------
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 
% Date last modified:
%
%
function h = qsmhub_setAllCallbacks(h)
set(h.TabGroup,                 'SelectionChangedFcn', {@test_Callback})
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
set(h.OneStop_pushbutton_start,	'Callback',{@PushbuttonOneStopStart_Callback});
set(h.checkbox_excludeMask,     'Callback',{@CheckboxExcludeMask_Callback});
set(h.checkbox_MEDI_smv,        'Callback',{@CheckboxMEDISMV_Callback});
set(h.checkbox_MEDI_lambda_csf, 'Callback',{@CheckboxMEDILambdaCSF_Callback});
set(h.edit_excludeMask,         'Callback',{@EditRange01_Callback});
set(h.edit_LBV_depth,           'Callback',{@EditMin_Callback});
set(h.edit_LBV_peel,            'Callback',{@EditNonNegative_Callback});
set(h.edit_LBV_tol,             'Callback',{@EditNonNegative_Callback});
set(h.edit_PDF_maxIter,        	'Callback',{@EditNonNegative_Callback});
set(h.edit_PDF_tol,             'Callback',{@EditNonNegative_Callback});
set(h.edit_PDF_padSize,      	'Callback',{@EditNonNegative_Callback});
set(h.edit_RESHARP_lambda,     	'Callback',{@EditNonNegative_Callback});
set(h.edit_RESHARP_radius,     	'Callback',{@EditNonNegative_Callback});
set(h.edit_SHARP_radius,     	'Callback',{@EditNonNegative_Callback});
set(h.edit_SHARP_threshold,    	'Callback',{@EditNonNegative_Callback});
set(h.edit_VSHARP_minRadius,   	'Callback',{@EditNonNegative_Callback});
set(h.edit_VSHARP_maxRadius,   	'Callback',{@EditVSHARPMaxRadius_Callback});
set(h.edit_VSHARPSTI_smvSize,   'Callback',{@EditNonNegative_Callback});
set(h.edit_iHARPERELLA_maxIter,	'Callback',{@EditNonNegative_Callback});
set(h.edit_TKD_threshold,       'Callback',{@EditRange01_Callback});
set(h.edit_cfs_lambda,          'Callback',{@EditNonNegative_Callback});
set(h.edit_cfs_lambda,          'Callback',{@EditNonNegative_Callback});
set(h.edit_Star_padSize,       	'Callback',{@EditNonNegative_Callback});
set(h.edit_iLSQR_lambda,       	'Callback',{@EditNonNegative_Callback});
set(h.edit_iLSQR_maxIter,      	'Callback',{@EditNonNegative_Callback});
set(h.edit_iLSQR_tol,           'Callback',{@EditNonNegative_Callback});
set(h.edit_STIiLSQR_maxIter,   	'Callback',{@EditNonNegative_Callback});
set(h.edit_STIiLSQR_padSize,   	'Callback',{@EditNonNegative_Callback});
set(h.edit_STIiLSQR_threshold, 	'Callback',{@EditNonNegative_Callback});
set(h.edit_STIiLSQR_tol1,   	'Callback',{@EditNonNegative_Callback});
set(h.edit_STIiLSQR_tol2,   	'Callback',{@EditNonNegative_Callback});
set(h.edit_FANSI_lambda,        'Callback',{@EditNonNegative_Callback});
set(h.edit_FANSI_maxIter,       'Callback',{@EditNonNegative_Callback});
set(h.edit_FANSI_mu,            'Callback',{@EditNonNegative_Callback});
set(h.edit_FANSI_tol,           'Callback',{@EditNonNegative_Callback});
end