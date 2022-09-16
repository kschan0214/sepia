%% Wrapper_SWISMWI_SMWI(headerAndExtraData,output_structure,algorParam)
%
% Input
% --------------
% headerAndExtraData : structure contains extra header info/data for the algorithm
% output_structure   : structure contains output file path and output NIfTI header
% algorParam         : structure contains fields with algorithm-specific parameter(s)
%
% Output
% --------------
%
% Description: This is a wrapper function to access FANSI for SEPIA
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 3 August 2022 (v1.1)
% Date modified: 
%
%
function Wrapper_SWISMWI_SMWI(headerAndExtraData,output_structure,algorParam)
sepia_universal_variables;

% get algorithm parameters
algorParam = check_and_set_algorithm_default(algorParam);
thres           = algorParam.swismwi.threshold;
m               = algorParam.swismwi.m;
ismIP           = algorParam.swismwi.ismIP;
slice_mIP       = algorParam.swismwi.slice_mIP;
isParamagnetic  = algorParam.swismwi.isParamagnetic;
isDiamagnetic   = algorParam.swismwi.isDiamagnetic;

% get extra data such as magnitude/weights/B0 direction/TE/etc.
headerAndExtraData                      = check_and_set_SEPIA_header_data(headerAndExtraData);
% b0dir = headerAndExtraData.sepia_header.B0_dir;
% matrixSize    = headerAndExtraData.sepia_header.matrixSize;
chi     = get_variable_from_headerAndExtraData(headerAndExtraData, 'phasechi');
magn    = get_variable_from_headerAndExtraData(headerAndExtraData, 'magnitude');

outputNiftiTemplate = output_structure.outputNiftiTemplate;
outputFileFullPrefix    = output_structure.outputFileFullPrefix;

%% SWI
disp('Computing SMWI...');

[pSMWI, dSMWI] = smwi(magn,chi,thres,m);

% minimum intensity projection
if ismIP
    psmwi_mIP = zeros(size(pSMWI));
    dsmwi_mIP = zeros(size(pSMWI));
    for kz = 1:size(pSMWI,3) - slice_mIP
        psmwi_mIP(:,:,kz,:) = min(pSMWI(:,:,kz:kz+slice_mIP-1,:),[],3);
        dsmwi_mIP(:,:,kz,:) = min(dSMWI(:,:,kz:kz+slice_mIP-1,:),[],3);
    end
end

%% save nii
disp('Saving SMWI results...');
if isParamagnetic
    save_nii_quick(outputNiftiTemplate,pSMWI, [outputFileFullPrefix 'smwi-paramagnetic.nii.gz']);
    if ismIP
        save_nii_quick(outputNiftiTemplate,psmwi_mIP, [outputFileFullPrefix 'smwi-mIP-paramagnetic.nii.gz']);
    end
end
if isDiamagnetic
    save_nii_quick(outputNiftiTemplate,dSMWI, [outputFileFullPrefix 'smwi-diamagnetic.nii.gz']);
    if ismIP
        save_nii_quick(outputNiftiTemplate,dsmwi_mIP, [outputFileFullPrefix 'smwi-mIP-diamagnetic.nii.gz']);
    end
end

disp('Done!');
end

%% set default parameter if not specified
function algorParam2 = check_and_set_algorithm_default(algorParam)

algorParam2 = algorParam;

try algorParam2.swismwi.threshold       = algorParam.swismwi.threshold;        catch; algorParam2.swismwi.threshold         = 1;	end
try algorParam2.swismwi.m            	= algorParam.swismwi.m;                catch; algorParam2.swismwi.m             	= 4;	end
try algorParam2.swismwi.ismIP         	= algorParam.swismwi.ismIP;            catch; algorParam2.swismwi.ismIP          	= true; end
try algorParam2.swismwi.slice_mIP       = algorParam.swismwi.slice_mIP;        catch; algorParam2.swismwi.slice_mIP      	= 4;	end
try algorParam2.swismwi.isParamagnetic	= algorParam.swismwi.isParamagnetic;   catch; algorParam2.swismwi.isParamagnetic	= true;	end
try algorParam2.swismwi.isDiamagnetic	= algorParam.swismwi.isDiamagnetic;    catch; algorParam2.swismwi.isDiamagnetic     = false;end

end
