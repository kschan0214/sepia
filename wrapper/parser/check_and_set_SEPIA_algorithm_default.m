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
% Date modified: 13 June 2020 (v0.8.0)
% Date modified: 12 September 2022 (v1.1.0)
%
%
function algorParam2 = check_and_set_SEPIA_algorithm_default(algorParam)

sepia_universal_variables;

algorParam2 = algorParam;

try algorParam2.general.isInvert   	= algorParam.general.isInvert; 	catch; algorParam2.general.isInvert	= false;	end
try algorParam2.general.isGPU       = algorParam.general.isGPU;     catch; algorParam2.general.isGPU    = false;    end
try algorParam2.general.isBET   	= algorParam.general.isBET; 	catch; algorParam2.general.isBET   	= false;	end
try algorParam2.general.fractional_threshold   	= algorParam.general.fractional_threshold; 	catch; algorParam2.general.fractional_threshold   	= 0.5;	end
try algorParam2.general.gradient_threshold   	= algorParam.general.gradient_threshold; 	catch; algorParam2.general.gradient_threshold   	= 0;	end

% default method is MEDI nonlinear fitting + Laplacian + no eddy correct + no voxel exclusion
try algorParam2.unwrap.echoCombMethod       = algorParam.unwrap.echoCombMethod;         catch; algorParam2.unwrap.echoCombMethod        = methodEchoCombineName{2};	end
% default phase unwrapping method is Laplacian
try algorParam2.unwrap.unwrapMethod         = algorParam.unwrap.unwrapMethod;           catch; algorParam2.unwrap.unwrapMethod          = methodUnwrapName{1};   	end
try algorParam2.unwrap.isEddyCorrect        = algorParam.unwrap.isEddyCorrect;          catch; algorParam2.unwrap.isEddyCorrect         = 0;                    end
try algorParam2.unwrap.excludeMaskThreshold	= algorParam.unwrap.excludeMaskThreshold;	catch; algorParam2.unwrap.excludeMaskThreshold	= Inf;                  end
% for the rest, if the parameter does not exist then initiates it with an empty array
try algorParam2.unwrap.isSaveUnwrappedEcho	= algorParam.unwrap.isSaveUnwrappedEcho;	catch; algorParam2.unwrap.isSaveUnwrappedEcho	= 0;                    end
try algorParam2.unwrap.excludeMethod        = algorParam.unwrap.excludeMethod;          catch; algorParam2.unwrap.excludeMethod         = 'Weighting map';      end
try algorParam2.unwrap.unit                 = algorParam.unwrap.unit;                   catch; algorParam2.unwrap.unit         = 'Hz';      end

% default background field removal method is VSHARP
try algorParam2.bfr.method              = algorParam.bfr.method;        catch; algorParam2.bfr.method = methodBFRName{5};   end
try algorParam2.bfr.erode_radius        = algorParam.bfr.erode_radius;	catch; algorParam2.bfr.erode_radius = 0;            end
try algorParam2.bfr.erode_before_radius	= algorParam.bfr.erode_before_radius;	catch; algorParam2.bfr.erode_before_radius = 0;            end
% try algorParam2.bfr.refine          = algorParam.bfr.refine;        catch; algorParam2.bfr.refine = false;              end
try algorParam2.bfr.refine_method       = algorParam.bfr.refine_method;	catch; algorParam2.bfr.refine_method = 'polyfit';	end
try algorParam2.bfr.refine_order        = algorParam.bfr.refine_order;	catch; algorParam2.bfr.refine_order = 4;         	end

% default background field removal method is TKD
try algorParam2.qsm.method              = algorParam.qsm.method;            catch; algorParam2.qsm.method           = methodQSMName{1};	end
try algorParam2.qsm.reference_tissue	= algorParam.qsm.reference_tissue;	catch; algorParam2.qsm.reference_tissue	= 'None';      	end

if strcmp(algorParam2.unwrap.echoCombMethod ,'ROMEO total field calculation')
    algorParam2.unwrap.unwrapMethod = 'ROMEO';
end

end