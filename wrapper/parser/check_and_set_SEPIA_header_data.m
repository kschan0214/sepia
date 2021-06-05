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
function headerAndExtraData2 = check_and_set_SEPIA_header_data(headerAndExtraData)

sepia_universal_variables;

headerAndExtraData2 = headerAndExtraData;

% try to get the data from input, if not then set a default value
% default: B0 align to z-direction
try headerAndExtraData2.b0dir       = headerAndExtraData.b0dir./norm(headerAndExtraData.b0dir);     catch; headerAndExtraData2.b0dir    = [0,0,1]; end 
% default 3 T
try headerAndExtraData2.b0          = headerAndExtraData.b0;        catch; headerAndExtraData2.b0       = 3; end
try headerAndExtraData2.te          = headerAndExtraData.te;        catch; headerAndExtraData2.te       = 40e-3; end
try headerAndExtraData2.delta_TE    = headerAndExtraData.delta_TE; 	catch; headerAndExtraData2.delta_TE = 40e-3; end
try headerAndExtraData2.CF          = headerAndExtraData.CF;        catch; headerAndExtraData2.CF       = headerAndExtraData2.b0*gyro*1e6; end

try headerAndExtraData2.weights     = headerAndExtraData.weights;   catch; headerAndExtraData2.weights   = []; end
try headerAndExtraData2.magn        = headerAndExtraData.magn;      catch; headerAndExtraData2.magn      = []; end
try headerAndExtraData2.phase       = headerAndExtraData.phase;     catch; headerAndExtraData2.phase     = []; end
try headerAndExtraData2.mask_ref    = headerAndExtraData.mask_ref; 	catch; headerAndExtraData2.mask_ref  = []; end
try headerAndExtraData2.initGuess   = headerAndExtraData.initGuess;	catch; headerAndExtraData2.initGuess = []; end
try headerAndExtraData2.N_std       = headerAndExtraData.N_std;     catch; headerAndExtraData2.N_std     = []; end


end