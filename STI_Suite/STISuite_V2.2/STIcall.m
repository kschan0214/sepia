params2.H=[0 0 1];
params.voxelsize= MyPARAMS.params.res;
params.niter=100;
params.TE=1;
params.B0=3;
params.tol_step1=0.01;
params.tol_step2=0.001;
params.Kthreshold=0.25;
params.padsize=4;
clear Confidence
clear M0
clear T2map
clear datarecon

SF=ScalingFactor(params.B0,params.TE);

[Susceptibility] = QSM_iLSQR(TissuePhase,NewMask,'params',params);

