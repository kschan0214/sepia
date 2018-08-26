function r2 = arlo(te,y)
% r2 = arlo(te,y)
%    
% output
%   r2 r2-map (in Hz)
%      empty when only one echo is provided
%
% input
%   te  array containing te values (in s)
%   y   a multi-echo data set of arbitrary dimension
%       echo should be the last dimension
%       
% If you use this function please cite
% 
% Pei M, Nguyen TD, Thimmappa ND, Salustri C, Dong F, Cooper MA, Li J, 
% Prince MR, Wang Y. Algorithm for fast monoexponential fitting based 
% on Auto-Regression on Linear Operations (ARLO) of data. 
% Magn Reson Med. 2015 Feb;73(2):843-50. doi: 10.1002/mrm.25137. 
% Epub 2014 Mar 24. PubMed PMID: 24664497; 
% PubMed Central PMCID:PMC4175304.

nte = length(te);
if nte<2
    r2 =[];
    return
end

sz=size(y);
edx = numel(sz);
if sz(edx)~=nte
    error(['Last dimension of y has size ' num2str(sz(edx)) ...
        ', expected ' num2str(nte) ]);
end

     yy=zeros(sz(1:end-1));
     yx=zeros(sz(1:end-1));
beta_yx=zeros(sz(1:end-1));
beta_xx=zeros(sz(1:end-1));
s1=[]; d1=[];
crd=repmat({':'}, [1 edx]);
crd0=crd;crd1=crd;crd2=crd;
for j=1:nte-2
    alpha = (te(j+2)-te(j))*(te(j+2)-te(j))/2/(te(j+1)-te(j));
    tmp = (2*te(j+2)*te(j+2) - te(j)*te(j+2) - te(j)*te(j) + 3*te(j)*te(j+1) -3*te(j+1)*te(j+2))/6; 
    beta = tmp/(te(j+2)-te(j+1)); 
    gamma = tmp/(te(j+1)-te(j));
    crd0{edx}=j;crd1{edx}=j+1;crd2{edx}=j+2;
%     [te(j+2)-te(j)-alpha+gamma alpha-beta-gamma beta]/((te(2)-te(1))/3)
    y1 = y(crd0{:})*(te(j+2)-te(j)-alpha+gamma)+y(crd1{:})*(alpha-beta-gamma)+y(crd2{:})*beta;
    x1 = y(crd0{:})-y(crd2{:});
         yy=yy+y1.*y1;
         yx=yx+y1.*x1;
    beta_yx=beta_yx+beta*y1.*x1;
    beta_xx=beta_xx+beta*x1.*x1;

end
r2 = (yx + beta_xx)./(beta_yx + yy);
r2(isnan(r2)) = 0;
r2(isinf(r2)) = 0;
