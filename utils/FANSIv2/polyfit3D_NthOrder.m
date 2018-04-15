function I_fitted = polyfit3D_NthOrder(I,mask,Order)
%polyfit the part of the image 'I' that is inside the mask

tic
matrix(:, 1) = ones(sum(mask(:)), 1);
FullMatrix(:, 1) = double(ones(size(mask,1)*size(mask,2)*size(mask,3),1));

x = (0 : size(mask,2) - 1)/size(mask,2) - 1/2;
Xmat = repmat(x, [size(mask, 1), 1, size(mask,3)]);
y = ((0 : size(mask,1) - 1)/size(mask,1) - 1/2).';
Ymat = repmat(y, [1, size(mask, 2),size(mask,3)]);
z = permute(((0 : size(mask,3) - 1)/size(mask,3) - 1/2),[1 3 2]);
Zmat = repmat(z, [size(mask,1), size(mask, 2),1]);

for CurrentOrder = 1:Order
    for Xpower = 0:CurrentOrder
        for Zpower = 0:(CurrentOrder-Xpower)
            matrix = [matrix, Xmat(mask == 1).^Xpower .*Zmat(mask == 1).^(Zpower) .*Ymat(mask == 1).^(CurrentOrder-(Xpower+Zpower)) ];
            FullMatrix = [FullMatrix, Xmat(:).^Xpower .*Zmat(:).^(Zpower) .*Ymat(:).^(CurrentOrder-(Xpower+Zpower))];
        end
    end
end

toc

I_fitted = zeros(size(I));
for c = 1:size(I,4)
    if mod(c,5)==0
        disp(num2str(c))
    end
    temp = I(:,:,:,c);
    Ivec = temp(mask==1);

    coefs = matrix \ Ivec; % determine coefficients based on data in W2
    I_fittedVec = FullMatrix * coefs;
    I_fitted(:,:,:,c) = reshape(I_fittedVec,size(I,1),size(I,2),size(I,3));
end
    
    
end
