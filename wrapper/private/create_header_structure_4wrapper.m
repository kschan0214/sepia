%% Standard script to save SEPIA header to the a variable structure for IOwrapper
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 22 Jan 2021
% Date modified: 
%
%
headerAndExtraData              = struct();
headerAndExtraData.b0           = B0;
headerAndExtraData.b0dir        = B0_dir;
headerAndExtraData.te           = TE;
headerAndExtraData.delta_TE     = delta_TE;
headerAndExtraData.CF           = CF;
headerAndExtraData.voxelSize    = voxelSize;
headerAndExtraData.matrixSize   = matrixSize;
