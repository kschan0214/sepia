% Projection onto Dipole Fields (PDF)
%   [p1, dp1, relres, p0]=Fit_ppm_complex_TE(M,TE)
%    
%   output
%   p1 - field map, may need further unwrapping
%   dp1 - a priori error estimate
%   relres - relative residual
%   p0 - initla phase
%
%   input
%   M - a multi-echo and could be a multi-channel dataset
%       echo needs to be the 4th dimension
%       channel needs to be the 5th dimension
%   TE - echo time in seconds
%
%   When using the code, please cite 
%   T. Liu et al. MRM 2013;69(2):467-76
%   B. Kressler et al. IEEE TMI 2010;29(2):273-81
%   de Rochefort et al. MRM 2008;60(4):1003-1009
%
%   The coil combination method is similar to
%   MA. Bernstein et al. MRM 1994;32:330-334
%
%   Adapted from a linear fitting created by Ludovic de Rochefort
%   Modified by Tian Liu on 2011.06.01
%   Last modified by Tian Liu on 2013.07.23


function [p1, dp1, relres, p0, iter]=Fit_ppm_complex_TE(M,TE)

if size(M,5)>1
% combine multiple coils together, assuming the coil is the fifth dimension
    M = sum(M.*conj( repmat(M(:,:,:,1,:),[1 1 1 size(M,4) 1])),5);  
    M = sqrt(abs(M)).*exp(1i*angle(M));
end

M= conj(M);
s0=size(M);
L_s0=length(s0);
nechos=size(M,L_s0);

M=reshape(M,[prod(s0(1:L_s0-1)),s0(L_s0)]);
s=size(M);

Y=angle(M(:,1:min(3,nechos)));
% estimate the slope
c=((Y(:,2)-Y(:,1)));
[m ind]=min([abs(c-2*pi),abs(c),abs(c+2*pi)],[],2);
c(ind==1)=c(ind==1)-2*pi;
c(ind==3)=c(ind==3)+2*pi;

% unwrap the second echo
cd=Y(:,2)-Y(:,1)-c;
Y(cd<-pi,2)=Y(cd<-pi,2)+2*pi;
Y(cd>pi,2)=Y(cd>pi,2)-2*pi;

% unwrap the third echo
cd=(Y(:,3)-Y(:,2))-(TE(3)-TE(2))/(TE(2)-TE(1))*c;
cd_minus=(cd<-pi);
Y(:,3)=Y(:,3)+cd_minus.*abs(fix((cd-pi)./(2*pi)))*2*pi;
cd_plus=(cd>pi);
Y(:,3)=Y(:,3)-cd_plus.*fix((cd+pi)./(2*pi))*2*pi;




A = [1  TE(1) ;1 TE(2);1 TE(3) ];
ip = A\Y(:,1:3)';
p0 = ip(1,:)';
p1 = ip(2,:)';

dp1 = p1;
tol = norm(p1(:))*1e-4;
iter = 0;
max_iter = 30;

% weigthed least square
% calculation of WA'*WA
v1=ones(1,nechos);
v2=reshape(TE,size(v1));
a11=sum(abs(M).^2.*(ones(s(1),1)*(v1.^2)),2);
a12=sum(abs(M).^2.*(ones(s(1),1)*(v1.*v2)),2);
a22=sum(abs(M).^2.*(ones(s(1),1)*(v2.^2)),2);
% inversion
d=a11.*a22-a12.^2;
ai11=a22./d;
ai12=-a12./d;
ai22=a11./d;

while ((norm(dp1)>tol) &&(iter<max_iter))
    iter = iter+1;
    W = abs(M).*exp(1i*(p0*v1 + p1*v2) );

    % projection
    pr1=sum(conj(1i*W).*(ones(s(1),1)*v1).*(M-W),2);
    pr2=sum(conj(1i*W).*(ones(s(1),1)*v2).*(M-W),2);

    dp0=real(ai11.*pr1+ai12.*pr2);
    dp1=real(ai12.*pr1+ai22.*pr2);
    dp1(isnan(dp1))=0;
    dp0(isnan(dp0))=0;
    
    %update
    p1 = p1+dp1;
    p0 = p0+dp0;
    

end

% error propagation
dp1=sqrt(ai22);
dp1(isnan(dp1)) = 0;
dp1(isinf(dp1)) = 0;

% relative residual
res = M - abs(M).*exp(1i*(p0*v1 + p1*v2) );
relres = sum(abs(res).^2,2)./sum(abs(M).^2,2);
relres(isnan(relres)) = 0;


p0=reshape(p0,s0(1:L_s0-1))*(TE(2)-TE(1));
p1=reshape(p1,s0(1:L_s0-1))*(TE(2)-TE(1));
dp1=reshape(dp1,s0(1:L_s0-1))*(TE(2)-TE(1));
relres = reshape(relres,s0(1:L_s0-1))*(TE(2)-TE(1));
p1(p1>pi)=mod(p1(p1>pi)+pi,2*pi)-pi;
p1(p1<-pi)=mod(p1(p1<-pi)+pi,2*pi)-pi;
    

