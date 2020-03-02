%% function algorParam2 = check_and_set_SEPIA_algorithm_default(algorParam)
%
% Input
% --------------
% algorParam    : strcuture variable contains algorithm parameters 
%
% Output
% --------------
% algorParam2   : strcuture variable contains all essential algorithm parameters 
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 28 Feb 2020
% Date last modified:
%
%
function algorParam2 = check_and_set_SEPIA_algorithm_default(algorParam)

algorParam2 = algorParam;

try algorParam2.general.isInvert   	= algorParam.general.isInvert; 	catch; algorParam2.general.isInvert	= false;	end
try algorParam2.general.isGPU       = algorParam.general.isGPU;     catch; algorParam2.general.isGPU    = false;    end
try algorParam2.general.isBET   	= algorParam.general.isBET; 	catch; algorParam2.general.isBET   	= false;	end
try algorParam2.general.fractional_threshold   	= algorParam.general.fractional_threshold; 	catch; algorParam2.general.fractional_threshold   	= 0.5;	end
try algorParam2.general.gradient_threshold   	= algorParam.general.gradient_threshold; 	catch; algorParam2.general.gradient_threshold   	= 0;	end

% default method is MEDI nonlinear fitting + Laplacian + no eddy correct + no voxel exclusion
try algorParam2.unwrap.echoCombMethod       = algorParam.unwrap.echoCombMethod;         catch; algorParam2.unwrap.echoCombMethod        = 'MEDI nonlinear fit';	end
% default phase unwrapping method is Laplacian
try algorParam2.unwrap.unwrapMethod         = algorParam.unwrap.unwrapMethod;           catch; algorParam2.unwrap.unwrapMethod          = 'Laplacian';          end
try algorParam2.unwrap.isEddyCorrect        = algorParam.unwrap.isEddyCorrect;          catch; algorParam2.unwrap.isEddyCorrect         = 0;                    end
try algorParam2.unwrap.excludeMaskThreshold	= algorParam.unwrap.excludeMaskThreshold;	catch; algorParam2.unwrap.excludeMaskThreshold	= Inf;                  end
% for the rest, if the parameter does not exist then initiates it with an empty array
try algorParam2.unwrap.subsampling          = algorParam.unwrap.subsampling;            catch; algorParam2.unwrap.subsampling           = [];                   end
try algorParam2.unwrap.isSaveUnwrappedEcho	= algorParam.unwrap.isSaveUnwrappedEcho;	catch; algorParam2.unwrap.isSaveUnwrappedEcho	= 0;                    end
try algorParam2.unwrap.excludeMethod        = algorParam.unwrap.excludeMethod;          catch; algorParam2.unwrap.excludeMethod         = 'Weighting map';      end

% default background field removal method is VSHARP
try algorParam2.bfr.method      = algorParam.bfr.method;        catch; algorParam2.bfr.method = 'vsharpsti';end
try algorParam2.bfr.radius      = algorParam.bfr.radius;        catch; algorParam2.bfr.radius = 10;         end
try algorParam2.bfr.refine      = algorParam.bfr.refine;        catch; algorParam2.bfr.refine = false;      end
try algorParam2.bfr.erode_radius= algorParam.bfr.erode_radius;	catch; algorParam2.bfr.erode_radius = 1;    end
% for the rest, if the parameter does not exist then initiates it with an empty array
try algorParam2.bfr.tol         = algorParam.bfr.tol;     	catch; algorParam2.bfr.tol = [];            end
try algorParam2.bfr.depth       = algorParam.bfr.depth;  	catch; algorParam2.bfr.depth = [];          end
try algorParam2.bfr.peel        = algorParam.bfr.peel;   	catch; algorParam2.bfr.peel = [];           end
try algorParam2.bfr.iteration   = algorParam.bfr.iteration;	catch; algorParam2.bfr.iteration = [];      end
try algorParam2.bfr.padSize     = algorParam.bfr.padSize; 	catch; algorParam2.bfr.padSize = [];        end
try algorParam2.bfr.alpha       = algorParam.bfr.alpha;    	catch; algorParam2.bfr.alpha = [];          end
try algorParam2.bfr.threshold   = algorParam.bfr.threshold;	catch; algorParam2.bfr.threshold = [];      end

% default background field removal method is TKD
try algorParam2.qsm.method      = algorParam.qsm.method;        catch; algorParam2.qsm.method       = 'tkd';	end
try algorParam2.qsm.threshold   = algorParam.qsm.threshold;     catch; algorParam2.qsm.threshold    = 0.15;     end
% for the rest, if the parameter does not exist then initiates it with an empty array
try algorParam2.qsm.radius      = algorParam.qsm.radius;        catch; algorParam2.qsm.radius       = [];       end
try algorParam2.qsm.lambda      = algorParam.qsm.lambda;        catch; algorParam2.qsm.lambda       = [];       end
try algorParam2.qsm.optimise   	= algorParam.qsm.optimise;     	catch; algorParam2.qsm.optimise     = [];       end
try algorParam2.qsm.maxiter   	= algorParam.qsm.maxiter;       catch; algorParam2.qsm.maxiter      = [];     	end
try algorParam2.qsm.tol1        = algorParam.qsm.tol1;          catch; algorParam2.qsm.tol1         = [];   	end
try algorParam2.qsm.tol2        = algorParam.qsm.tol2;          catch; algorParam2.qsm.tol2         = [];      	end
try algorParam2.qsm.padsize     = algorParam.qsm.padsize;       catch; algorParam2.qsm.padsize      = [];   	end
try algorParam2.qsm.tol         = algorParam.qsm.tol;           catch; algorParam2.qsm.tol          = [];      	end
try algorParam2.qsm.mu1         = algorParam.qsm.mu1;           catch; algorParam2.qsm.mu1          = [];      	end
try algorParam2.qsm.mu2         = algorParam.qsm.mu2;           catch; algorParam2.qsm.mu2          = [];    	end
try algorParam2.qsm.solver    	= algorParam.qsm.solver;    	catch; algorParam2.qsm.solver       = [];    	end
try algorParam2.qsm.constraint 	= algorParam.qsm.constraint;   	catch; algorParam2.qsm.constraint   = [];       end
try algorParam2.qsm.gradient_mode	= algorParam.qsm.gradient_mode;	catch; algorParam2.qsm.gradient_mode    = [];      	end
try algorParam2.qsm.isWeakHarmonic  = algorParam.qsm.isWeakHarmonic;catch; algorParam2.qsm.isWeakHarmonic	= [];      	end
try algorParam2.qsm.beta            = algorParam.qsm.beta;          catch; algorParam2.qsm.beta             = [];      	end
try algorParam2.qsm.muh             = algorParam.qsm.muh;           catch; algorParam2.qsm.muh              = [];      	end
try algorParam2.qsm.wData       = algorParam.qsm.wData;         catch; algorParam2.qsm.wData        = [];       end
try algorParam2.qsm.wGradient  	= algorParam.qsm.wGradient;    	catch; algorParam2.qsm.wGradient    = [];       end
try algorParam2.qsm.zeropad  	= algorParam.qsm.zeropad;    	catch; algorParam2.qsm.zeropad      = [];      	end
try algorParam2.qsm.isSMV    	= algorParam.qsm.isSMV;         catch; algorParam2.qsm.isSMV        = [];    	end
try algorParam2.qsm.isLambdaCSF	= algorParam.qsm.isLambdaCSF;	catch; algorParam2.qsm.isLambdaCSF  = [];   	end
try algorParam2.qsm.lambdaCSF 	= algorParam.qsm.lambdaCSF;    	catch; algorParam2.qsm.lambdaCSF    = [];     	end
try algorParam2.qsm.merit       = algorParam.qsm.merit;         catch; algorParam2.qsm.merit        = [];      	end
try algorParam2.qsm.stepSize   	= algorParam.qsm.stepSize;    	catch; algorParam2.qsm.stepSize    	= [];      	end
try algorParam2.qsm.percentage  = algorParam.qsm.percentage;    catch; algorParam2.qsm.percentage   = [];      	end

end