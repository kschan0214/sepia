%% Standard script to display input sepia header information
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 22 Jan 2021
% Date modified: 
%
%
disp('----------------------');
disp('Basic data information');
disp('----------------------');
fprintf('Voxel size(x,y,z)   = %s mm x %s mm x %s mm\n' ,num2str(voxelSize(1)),num2str(voxelSize(2)),num2str(voxelSize(3)));
fprintf('Matrix size(x,y,z)  = %s x %s x %s\n'          ,num2str(matrixSize(1)),num2str(matrixSize(2)),num2str(matrixSize(3)));
fprintf('B0 direction(x,y,z) = [%s; %s; %s]\n'          ,num2str(B0_dir(1)),num2str(B0_dir(2)),num2str(B0_dir(3)));
fprintf('Field strength      = %s T\n'                  ,num2str(B0));
fprintf('Number of echoes    = %s\n'                    ,num2str(length(TE)));
fprintf('TE1/dTE             = %s/%s ms\n'              ,num2str(TE(1)*1e3),num2str(delta_TE*1e3));
