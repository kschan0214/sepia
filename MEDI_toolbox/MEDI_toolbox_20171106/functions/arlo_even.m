function r2 = arlo(te,y)
% Pei M, Nguyen TD, Thimmappa ND, Salustri C, Dong F, Cooper MA, Li J, 
% Prince MR, Wang Y. Algorithm for fast monoexponential fitting based 
% on Auto-Regression on Linear Operations (ARLO) of data. 
% Magn Reson Med. 2015 Feb;73(2):843-50. doi: 10.1002/mrm.25137. 
% Epub 2014 Mar 24. PubMed PMID: 24664497; 
% PubMed Central PMCID:PMC4175304.

nte = length(te);
esp = te(2)-te(1);
sz=size(y);
edx = numel(sz);
if sz(edx)~=nte
    error(['Last dimension of y has size ' num2str(sz(edx)) ...
        ', expected ' num2str(nte) ]);
end
x1 = zeros([sz(1:end-1) nte-2]);
y1 = x1;

crd=repmat({':'}, [1 edx]);
crd0=crd;crd1=crd;crd2=crd;
crd0{edx}=1:nte-2;crd1{edx}=2:nte-1;crd2{edx}=3:nte;
y1=esp/3*(y(crd0{:})+4*y(crd1{:})+y(crd2{:}));
x1=y(crd0{:})-y(crd2{:});
r2 = (esp/3*sum(x1.^2,edx) + sum(y1.*x1,edx))./ ...
    (sum(y1.^2,edx) + esp/3*sum(y1.*x1,edx));
r2(isnan(r2)) = 0;
r2(isinf(r2)) = 0;
