%% sepiatest.compare_to_reference(testCase, actual, reference, tol, label)
%
% Description: compares a sepiatest.mask_stats() struct against a
% previously-saved reference struct (same shape), using testCase.verifyEqual
% with a relative+absolute tolerance, and a normalized-RMSE check on the
% coarse downsampled array. Uses verifyXXX (not assertXXX) so a single
% failing statistic is reported clearly without aborting the rest of the
% comparison (e.g. mean off but median fine is informative).
%
% Input
% --------------
% testCase  : a matlab.unittest.qualifications.Qualifiable (i.e. a TestCase)
% actual    : struct from sepiatest.mask_stats(...) computed on the fresh run
% reference : struct with the same fields, loaded from a references/*.mat file
% tol       : struct with fields relTol, absTol (see sepiatest.tolerance_for_method)
% label     : short string identifying this comparison in failure messages
%             (e.g. 'QSM:TKD')
%
function compare_to_reference(testCase, actual, reference, tol, label)

if nargin < 5 || isempty(label)
    label = 'comparison';
end

scalarFields = {'mean','std','median','p1','p5','p25','p75','p95','p99'};
for k = 1:numel(scalarFields)
    f = scalarFields{k};
    testCase.verifyEqual(actual.(f), reference.(f), ...
        'RelTol', tol.relTol, 'AbsTol', tol.absTol, ...
        sprintf('%s: stat "%s" outside tolerance (actual=%.6g, reference=%.6g, relTol=%.4g, absTol=%.4g)', ...
                 label, f, actual.(f), reference.(f), tol.relTol, tol.absTol));
end

testCase.verifyEqual(actual.nNaN, reference.nNaN, ...
    sprintf('%s: number of NaN voxels inside mask changed (actual=%d, reference=%d)', label, actual.nNaN, reference.nNaN));
testCase.verifyEqual(actual.nInf, reference.nInf, ...
    sprintf('%s: number of Inf voxels inside mask changed (actual=%d, reference=%d)', label, actual.nInf, reference.nInf));

% coarse spatial-pattern check via normalized RMSE on the downsampled array
refDs = reference.downsampled;
actDs = actual.downsampled;
if isequal(size(refDs), size(actDs))
    denom = max(abs(refDs(:)));
    if denom == 0
        denom = 1;
    end
    nrmse = sqrt(mean((actDs(:) - refDs(:)).^2)) / denom;
    testCase.verifyLessThanOrEqual(nrmse, max(tol.relTol, 0.05), ...
        sprintf('%s: coarse spatial-pattern NRMSE %.4g exceeds tolerance %.4g', label, nrmse, max(tol.relTol,0.05)));
else
    testCase.verifyFail(sprintf('%s: downsampled array size changed (actual=[%s], reference=[%s])', ...
        label, num2str(size(actDs)), num2str(size(refDs))));
end

end
