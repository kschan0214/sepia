%% tf = qsm_method_uses_mask_in_inversion(qsmMethod, solver)
%
% Input
% --------------
% qsmMethod : QSM dipole-inversion method name (algorParam.qsm.method /
%             methodQSMName entry)
% solver    : (optional) method-specific solver name, for methods with
%             multiple solvers (e.g. the "MRI Suscep. Calc." addon's
%             algorParam.qsm.solver)
%
% Output
% --------------
% tf : false for QSM dipole-inversion configurations known to apply the
%      mask only as a final multiplication *after* a closed-form/direct
%      k-space inversion (Truncated K-space Division, Direct Tikhonov) -
%      the mask never enters the deconvolution itself, so two-pass masking
%      has no effect on their output by construction. true otherwise (the
%      default assumption for iterative/regularised methods - e.g. FANSI,
%      MEDI, Iterative Tikhonov - where the mask enters the solve directly
%      and a refined mask can change the reconstructed values).
%
% Description: used to warn users in the GUI (and could be used to warn in
% batch/config-file runs) when they pair two-pass masking with a QSM
% method that cannot benefit from it. See
% sepia.documentation/docs/method/qsm/Two-pass-masking.rst for the
% method-dependence background and how this was established empirically.
%
% Kwok-shing Chan
% Date created: 10 September 2026
%
function tf = qsm_method_uses_mask_in_inversion(qsmMethod, solver)

if nargin < 2; solver = ''; end

tf = true;

if strcmpi(qsmMethod, 'TKD')
    tf = false;
elseif strcmpi(qsmMethod, 'MRI Suscep. Calc.') && ...
       any(strcmpi(solver, {'Truncated kspace division','Direct Tikhonov'}))
    tf = false;
end

end
