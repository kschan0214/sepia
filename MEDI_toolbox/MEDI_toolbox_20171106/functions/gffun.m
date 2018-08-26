function [g1, g2, in] = gffun(echo_laplacian,Mask)

matrix_size = size(Mask);
A = zeros(matrix_size);
B = zeros(matrix_size);

g1 = zeros([size(echo_laplacian,4), 1]);
g2 = zeros([size(echo_laplacian,4), 1]);
in = zeros([size(echo_laplacian,4), 1]);

for i = 1 : matrix_size(1)
    A(i,:,:) = i;
end
for i = 1 : matrix_size(2)
    B(:,i,:) = i;
end

for jj = 1 : size(echo_laplacian,4)
    
    indX = A(Mask>0);
    indY = B(Mask>0);
    V = squeeze(echo_laplacian(:,:,:,jj));
    V = V(Mask>0);
  
    X = [indX indY ones(length(indX),1)];
    b1 = X\V;
    
    g1(jj) = b1(1);
    g2(jj) = b1(2);
    in(jj) = b1(3);
end

end