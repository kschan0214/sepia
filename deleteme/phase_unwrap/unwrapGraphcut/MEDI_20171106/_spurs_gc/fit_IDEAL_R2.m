function [water fat freq R2s iter model]= fit_IDEAL_R2(s0, t, f_fat, f0, R2s, max_iter)
matrix_size = size(s0);
numvox = prod(matrix_size(1:end-1));
numte = matrix_size(end);

if nargin<6
    max_iter = 30;
end
if nargin <5
    R2s = zeros([1 numvox]);
end
if nargin<4;
    f0 = zeros([1 numvox]);
end
if nargin<3;
    f_fat = -420;
end

if numel(f0) == 1
    f0 = f0*ones([1 numvox]);
end
if numel(R2s) == 1
    R2s = R2s*ones([1 numvox]);
end

s0 = permute(s0, [length(matrix_size) 1:length(matrix_size)-1]);
s0 = reshape(s0,[numte numvox]);
R2s = reshape(R2s,[1 numvox]);
f0 = reshape(f0,[1 numvox]);
t = reshape(t,[numte 1]);
t = repmat(t,[1 numvox]);

% O = ones([numte numvox]).*exp(repmat(-R2s,[numte 1]).*t);
% C = real(exp(-1i*2*pi*f_fat*t)).*exp(repmat(-R2s, [numte 1]).*t);

O = ones([numte numvox]);
C = exp(-1i*2*pi*f_fat*t);


y = zeros([3 numvox]);
dy = zeros([3 numvox]);
dy(1,:) = 1e4;
iter = 0;

y(1,:) = f0+1i*R2s/(2*pi);

P = exp(-1i*2*pi*repmat(y(1,:),[numte 1]).*t);  %complex phasor
y(2:3,:) = invA(P.*O, P.*C, s0);

update = dy(1,:);
while (iter<max_iter)&&(sqrt(sum(real(update).^2,2)/numvox)>0.1)
    sn = P.*O.*repmat(y(2,:),[numte 1]) + P.*C.*repmat(y(3,:),[numte 1]);
    sr = s0 - sn;

%     gr = -2*pi*t.*(-repmat(y(2,:),[numte 1]).*z - repmat(y(3,:),[numte 1]).*O - repmat(y(4,:),[numte 1]).*d - repmat(y(5,:),[numte 1]).*C);
%     gi = -2*pi*t.*(repmat(y(2,:),[numte 1]).*O - repmat(y(3,:),[numte 1]).*z + repmat(y(4,:),[numte 1]).*C - repmat(y(5,:),[numte 1]).*d);

    Bcol01 = -1i*2*pi*t.*sn;
    dy = invB(Bcol01, P.*O, P.*C, sr);    
    y = y+dy;
    iter = iter+1;

    temp = y(1,:);
    temp(abs(imag(temp))*2*pi>50)=real(temp(abs(imag(temp))*2*pi>50));
%     y(1,:) = temp;  hann_low
    
    P = exp(-1i*2*pi*repmat(y(1,:),[numte 1]).*t);  %complex phasor
    y(2:3,:) = invA(P.*O, P.*C, s0);

    update = dy(1,:);
    update(isnan(update)) = 0;
    update(isinf(update)) = 0;
end

freq = reshape(real(y(1,:)),matrix_size(1:end-1));
R2s = -reshape(imag(y(1,:))*2*pi,matrix_size(1:end-1));
water = reshape(y(2,:),matrix_size(1:end-1));
fat = reshape(y(3,:),matrix_size(1:end-1));
model = P.*(O.*repmat(y(2,:),[numte 1]) + C.*repmat(y(3,:),[numte 1]));
model = reshape(model,matrix_size);
% s_model = [s0;model];
% figure; plot(s0,'ro'); hold on; plot(model,'bx');hold off;
% axis([-max(abs(real(s_model))) max(abs(real(s_model))) -max(abs(imag(s_model))) max(abs(imag(s_model)))]*1.2)

freq(isinf(freq)) = 0;
freq(isnan(freq)) = 0;
freq(abs(freq)>10e4) = 0;
R2s(isinf(R2s)) = 0;
R2s(isnan(R2s)) = 0;
water(isinf(water)) = 0;
water(isnan(water)) = 0;
water(abs(water)>10e5) = 0;
fat(isinf(fat)) = 0;
fat(isnan(fat)) = 0;
fat(abs(fat)>10e5) = 0;
model(isinf(model)) = 0;
model(isnan(model)) = 0;

% freq(isinf(freq)) = 0;
% freq(isnan(freq)) = 0;
% 
% R2s(isinf(R2s)) = 0;
% R2s(isnan(R2s)) = 0;
% water(isinf(water)) = 0;
% water(isnan(water)) = 0;
% 
% fat(isinf(fat)) = 0;
% fat(isnan(fat)) = 0;
% 
% model(isinf(model)) = 0;
% model(isnan(model)) = 0;


function x=invA(col1, col2, y)
% assemble A^H*A
a11 = sum(conj(col1).*col1, 1);
a12 = sum(conj(col1).*col2, 1);
a22 = sum(conj(col2).*col2, 1);

% inversion of A^H*A
d = (a11.*a22 - a12.*conj(a12));
ia11 = a22./d;
ia12 = -a12./d;
ia22 = a11./d;

% y project onto A^H
py1 = sum(conj(col1).*y,1);
py2 = sum(conj(col2).*y,1);
% calculate x
x(1,:) = sum(ia11.*py1 + ia12.*py2, 1);
x(2,:) = sum(conj(ia12).*py1 + ia22.*py2, 1);

function x=invB(col1, col2, col3, y)
% assemble B^H*B
b11 = sum(conj(col1).*col1,1);
b12 = sum(conj(col1).*col2,1);
b13 = sum(conj(col1).*col3,1);
b22 = sum(conj(col2).*col2,1);
b23 = sum(conj(col2).*col3,1);
b33 = sum(conj(col3).*col3,1);

% inversion of B'*B
d = (b13.*conj(b12).*conj(b23) + b11.*b22.*b33 + b12.*b23.*conj(b13) - b13.*b22.*conj(b13) - b11.*b23.*conj(b23) - b12.*b33.*conj(b12));
ib11 = (b22.*b33 - b23.*conj(b23))./d;
ib12 = -(b12.*b33 - b13.*conj(b23))./d;
ib13 = (b12.*b23 - b13.*b22)./d;
ib22 = (b11.*b33 - b13.*conj(b13))./d;
ib23 = -(b11.*b23 - b13.*conj(b12))./d;
ib33 =  (b11.*b22 - b12.*conj(b12))./d;

% y project onto B'
py1 = sum(conj(col1).*y,1);
py2 = sum(conj(col2).*y,1);
py3 = sum(conj(col3).*y,1);
% calculate x
x(1,:) = sum(ib11.*py1 + ib12.*py2 + ib13.*py3, 1);
x(2,:) = sum(conj(ib12).*py1 + ib22.*py2 + ib23.*py3, 1);
x(3,:) = sum(conj(ib13).*py1 + conj(ib23).*py2 + ib33.*py3, 1);
