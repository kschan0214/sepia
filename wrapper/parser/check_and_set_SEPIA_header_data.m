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
% Date modified: 13 August 2021 (v1.0)
%
%
function headerAndExtraData2 = check_and_set_SEPIA_header_data(headerAndExtraData)

sepia_universal_variables;

headerAndExtraData2 = headerAndExtraData;

% try to get the data from input, if not then set a default value
% default: B0 align to z-direction
try headerAndExtraData2.sepia_header.B0_dir    	= headerAndExtraData.sepia_header.B0_dir./norm(headerAndExtraData.sepia_header.B0_dir);     catch; headerAndExtraData2.sepia_header.B0_dir    = [0,0,1]; end 
% default 3 T
try headerAndExtraData2.sepia_header.B0          = headerAndExtraData.sepia_header.B0;      catch; headerAndExtraData2.sepia_header.B0       = 3; end
try headerAndExtraData2.sepia_header.TE          = headerAndExtraData.sepia_header.TE;      catch; headerAndExtraData2.sepia_header.TE       = 40e-3; end
try headerAndExtraData2.sepia_header.delta_TE    = headerAndExtraData.sepia_header.delta_TE;catch; headerAndExtraData2.sepia_header.delta_TE = 40e-3; end
try headerAndExtraData2.sepia_header.CF          = headerAndExtraData.sepia_header.CF;   	catch; headerAndExtraData2.sepia_header.CF       = headerAndExtraData2.sepia_header.b0*gyro*1e6; end

try headerAndExtraData2.weights     = headerAndExtraData.weights;   catch; headerAndExtraData2.weights   = []; end
try headerAndExtraData2.magnitude 	= headerAndExtraData.magnitude;	catch; headerAndExtraData2.magnitude = []; end
try headerAndExtraData2.phase       = headerAndExtraData.phase;     catch; headerAndExtraData2.phase     = []; end
try headerAndExtraData2.mask_ref    = headerAndExtraData.mask_ref; 	catch; headerAndExtraData2.mask_ref  = []; end
try headerAndExtraData2.initGuess   = headerAndExtraData.initGuess;	catch; headerAndExtraData2.initGuess = []; end
try headerAndExtraData2.fieldmapSD 	= headerAndExtraData.fieldmapSD;     catch; headerAndExtraData2.fieldmapSD     = []; end

try headerAndExtraData2.availableFileList.phase = headerAndExtraData.availableFileList.phase;           catch; headerAndExtraData2.availableFileList.phase = []; end
try headerAndExtraData2.availableFileList.magnitude = headerAndExtraData.availableFileList.magnitude;   catch; headerAndExtraData2.availableFileList.magnitude = []; end
try headerAndExtraData2.availableFileList.mask = headerAndExtraData.availableFileList.mask;             catch; headerAndExtraData2.availableFileList.mask = []; end
try headerAndExtraData2.availableFileList.fieldmapSD = headerAndExtraData.availableFileList.fieldmapSD; catch; headerAndExtraData2.availableFileList.fieldmapSD = []; end
try headerAndExtraData2.availableFileList.weights = headerAndExtraData.availableFileList.weights; catch; headerAndExtraData2.availableFileList.weights = []; end

try headerAndExtraData2.phasechi 	= headerAndExtraData.phasechi;      catch; headerAndExtraData2.phasechi     = []; end
try headerAndExtraData2.availableFileList.phasechi 	= headerAndExtraData.availableFileList.phasechi;     catch; headerAndExtraData2.availableFileList.phasechi     = []; end

end