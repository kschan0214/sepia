%% Standard script to display input sepia header information
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 22 Jan 2021
% Date modified: 12 August 2021 (v1.0)
%
%
disp('----------------------');
disp('Basic data information');
disp('----------------------');
fprintf('Voxel size(x,y,z)   = %s mm x %s mm x %s mm\n' ,num2str(sepia_header.voxelSize(1)),num2str(sepia_header.voxelSize(2)),num2str(sepia_header.voxelSize(3)));
fprintf('Matrix size(x,y,z)  = %s x %s x %s\n'          ,num2str(sepia_header.matrixSize(1)),num2str(sepia_header.matrixSize(2)),num2str(sepia_header.matrixSize(3)));
fprintf('B0 direction(x,y,z) = [%s; %s; %s]\n'          ,num2str(sepia_header.B0_dir(1)),num2str(sepia_header.B0_dir(2)),num2str(sepia_header.B0_dir(3)));
fprintf('Field strength      = %s T\n'                  ,num2str(sepia_header.B0));
fprintf('Number of echoes    = %s\n'                    ,num2str(length(sepia_header.TE)));
fprintf('TE1/dTE             = %s/%s ms\n'              ,num2str(sepia_header.TE(1)*1e3),num2str(sepia_header.delta_TE*1e3));
