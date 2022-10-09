%% [chi] = Wrapper_QSM_MEDI(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
%
% Input
% --------------
% localField    : local field map (tissue fields), in Hz
% mask          : signal mask
% matrixSize    : size of the input image
% voxelSize     : spatial resolution of each dimension of the data, in mm
% algorParam    : structure contains fields with algorithm-specific parameter(s)
% headerAndExtraData : structure contains extra header info/data for the algorithm
%
% Output
% --------------
% chi           : magnetic susceptibility map, in ppm
%
% Description: This is a wrapper function to access MEDI for SEPIA
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 8 March 2020
% Date modified: 19 Jan 2021
% Date modified: 13 August 2021 (v1.0)
%
%
function [chi] = Wrapper_QSM_MEDI(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
sepia_universal_variables;

% get algorithm parameters
algorParam = check_and_set_algorithm_default(algorParam);
method      = algorParam.qsm.method;
lambda      = algorParam.qsm.lambda;
wData       = algorParam.qsm.wData;
pad         = algorParam.qsm.zeropad;
percentage  = algorParam.qsm.percentage;
isSMV       = algorParam.qsm.isSMV;
radius      = algorParam.qsm.radius;
isMerit     = algorParam.qsm.merit;
isLambdaCSF = algorParam.qsm.isLambdaCSF;
lam_CSF     = algorParam.qsm.lambdaCSF;

wGrad       = 1;

% get extra data such as magnitude/weights/B0 direction/TE/etc.
headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);
B0_dir      = headerAndExtraData.sepia_header.B0_dir;
delta_TE    = headerAndExtraData.sepia_header.delta_TE;
TE          = headerAndExtraData.sepia_header.TE;
CF          = headerAndExtraData.sepia_header.CF;
N_std       = get_variable_from_headerAndExtraData(headerAndExtraData, 'weights', matrixSize); % contrast of N_std is inverse of weights, still need processing
iMag        = get_variable_from_headerAndExtraData(headerAndExtraData, 'magnitude', matrixSize);

matrix_size = matrixSize;
voxel_size  = voxelSize;
RDF         = localField;
Mask        = mask;
iFreq       = [];
Mask_CSF    = [];

% add path
sepia_addpath('MEDI');
addpath(fullfile(SEPIA_HOME,'misc','qsm_algorithm','MEDI_L1'));

%% Preparation
N_std               = (1./N_std).*mask;
N_std(isnan(N_std)) = 0;
N_std(isinf(N_std)) = 0;
N_std = N_std / norm(N_std(mask>0));

if length(TE) < 3 && isLambdaCSF
    warning('CSF regularisation requires data with more than 3 echoes');
    disp('No CSF regularisation will be used.');
    isLambdaCSF = false;
end

% zero reference using CSF requires CSF mask
if isLambdaCSF && ~isempty(iMag)
    disp('Extracting CSF mask....');
    % R2* mapping
    r2s         = arlo(TE,iMag);
    Mask_CSF    = extract_CSF(r2s,mask,voxelSize)>0;
%             magn    = sqrt(sum(magn.^2,4));
end

% compress 4D info to 3D
iMag = sqrt(sum(abs(iMag).^2,4));
        
% MEDI input expects local field in rad
RDF = RDF*2*pi*delta_TE;

tmp_output_dir  = [pwd filesep];
tmp_filename    = [tmp_output_dir 'tmp_RDF.mat'];
% matching naming convension for MEDI_L1

save(tmp_filename,'iFreq','RDF','N_std','iMag','Mask','matrix_size','voxel_size','delta_TE','CF','B0_dir','Mask_CSF');

if exist(fullfile('.','results'),'dir')
    isResultDirMEDIExist = true;
else
    isResultDirMEDIExist = false;
end

%% Display algorithm parameters
disp('The following parameters are being used...');
disp(['Regularisation lambda	= ' num2str(lambda)]);
disp(['Data mode                = ' num2str(wData)]);
disp(['Zero padding (x,y,z)     = ' num2str(pad)]);
disp(['Mask precentage (%)      = ' num2str(percentage)]);
disp(['Is SMV?                  = ' num2str(isSMV)]);
if isSMV
    disp(['Radius                   = ' num2str(radius)]);
end
disp(['Is Merit?                = ' num2str(isMerit)]);
disp(['Use CSF regularisation?  = ' num2str(isLambdaCSF)]);
if isLambdaCSF
    disp(['CSF lambda               = ' num2str(lam_CSF)])
end

%% main
percentage = percentage/100;

if isSMV
    chi = MEDI_L1('filename',tmp_filename,'lambda',lambda,'data_weighting',wData,'gradient_weighting',wGrad,...
              'merit',isMerit,'smv',radius,'zeropad',pad,'lambda_CSF',lam_CSF,'percentage',percentage);
    SphereK = single(sphere_kernel(matrix_size, voxel_size,radius));
    mask = SMV(Mask, SphereK)>0.999;
else
    chi = MEDI_L1('filename',tmp_filename,'lambda',lambda,'data_weighting',wData,'gradient_weighting',wGrad,...
              'merit',isMerit,'zeropad',pad,'lambda_CSF',lam_CSF,'percentage',percentage);
end

chi = chi .* mask;

% clean up MEDI output and temp files 
delete(tmp_filename);
if isResultDirMEDIExist
    fileno=getnextfileno(['results' filesep],'x','.mat') - 1;
    resultsfile=strcat(['results' filesep 'x'],sprintf('%08u',fileno), '.mat');
    delete(resultsfile)
else
    rmdir(fullfile(pwd,'results'),'s');
end

end

%% set default parameter for unspecific input
function algorParam2 = check_and_set_algorithm_default(algorParam)

algorParam2 = algorParam;

try algorParam2.qsm.lambda      = algorParam.qsm.lambda;        catch; algorParam2.qsm.lambda       = 1000; end
try algorParam2.qsm.wData       = algorParam.qsm.wData;         catch; algorParam2.qsm.wData        = 1;    end
try algorParam2.qsm.zeropad     = algorParam.qsm.zeropad;       catch; algorParam2.qsm.zeropad      = 0;    end
try algorParam2.qsm.percentage  = algorParam.qsm.percentage;    catch; algorParam2.qsm.percentage   = 0.9;  end
try algorParam2.qsm.isSMV       = algorParam.qsm.isSMV;         catch; algorParam2.qsm.isSMV        = false;end
try algorParam2.qsm.radius      = algorParam.qsm.radius;        catch; algorParam2.qsm.radius       = 5;    end
try algorParam2.qsm.merit       = algorParam.qsm.merit;         catch; algorParam2.qsm.merit        = false;end
try algorParam2.qsm.isLambdaCSF = algorParam.qsm.isLambdaCSF;   catch; algorParam2.qsm.isLambdaCSF  = false;end
try algorParam2.qsm.lambdaCSF   = algorParam.qsm.lambdaCSF;     catch; algorParam2.qsm.lambdaCSF    = 100;  end

end