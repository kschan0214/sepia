function x_cr = sparsa_solver_kspace_SL0(x0,param,tau)

%%  work well for brain with weigthing in phase
% A = @(x) param.PhiW*(param.iFDF*x);
% AT =@(x) param.iFDF'*(param.PhiW*x);
%  y = param.PhiW*param.phase;
%%
A = @(x) param.iFDF*x;
AT = @(x) param.iFDF'*x;
y = param.phase;
phi= @(x) param.TV*x;

 sigma_min=0.004;
 sigma_decrease_factor=0.5;
 mu_0 =2;
 L=3;
 
x_cr = SL0(A, y, sigma_min, sigma_decrease_factor, mu_0, L, AT);





