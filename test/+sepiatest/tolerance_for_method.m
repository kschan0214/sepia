%% tol = sepiatest.tolerance_for_method(category, methodName)
%
% Description: returns the numeric tolerance to use when comparing a
% method's output stats (sepiatest.mask_stats) against its saved
% reference (sepiatest.compare_to_reference). Deterministic closed-form
% methods get a tight tolerance; iterative solvers (whose exact iteration
% count/convergence can vary slightly across MATLAB/BLAS/FFT library
% versions) get a looser one. Unrecognised methods fall back to a
% conservative default rather than erroring, so a newly-added method is
% still testable before someone remembers to categorise it here.
%
% Output
% --------------
% tol : struct with fields relTol (relative tolerance) and absTol
%       (absolute floor, used when the reference value is near zero)
%
function tol = tolerance_for_method(category, methodName)

switch lower(category)

    case 'qsm'
        tight = {'tkd','closed-form solution'};
        loose = {'medi','fansi','star-qsm','ilsqr','sti suite ilsqr','ndi','lsqr+heidi'};
        if any(strcmpi(methodName, tight))
            tol = struct('relTol',1e-4,'absTol',1e-6);
        elseif any(strcmpi(methodName, loose))
            tol = struct('relTol',0.05,'absTol',1e-4);
        else
            tol = struct('relTol',0.05,'absTol',1e-4); % conservative default (deep-learning/add-on methods, etc.)
        end

    case 'bfr'
        tol = struct('relTol',0.02,'absTol',1e-6);

    case 'unwrap'
        tol = struct('relTol',1e-3,'absTol',1e-6);

    case 'r2s'
        tol = struct('relTol',1e-3,'absTol',1e-3);

    case 'twopass'
        tol = struct('relTol',0.05,'absTol',1e-4);

    otherwise
        tol = struct('relTol',0.05,'absTol',1e-4);
end

end
