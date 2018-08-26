function phasor = phasorprep(slope1, slope2, intrcp, matrix_size)

A = zeros(matrix_size);
B = zeros(matrix_size);
for i = 1 : matrix_size(1)
    A(i,:,:) = i;
end
for i = 1 : matrix_size(2)
    B(:,i,:) = i;
end
phasor = zeros(matrix_size);
phasor = repmat(phasor, [1 1 1 length(slope1)]);

for i = 1 : size(phasor,4)
    phasor(:,:,:,i) = A*slope1(i) + B*slope2(i) + intrcp(i);
end

end