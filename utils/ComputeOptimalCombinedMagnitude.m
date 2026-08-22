function combinedMagnitude = ComputeOptimalCombinedMagnitude(te,r2s,magn)

% Compute weights
wTE = zeros(size(magn));
for kt=1:length(te)
wTE(:,:,:,kt) = te(kt)*exp(-r2s.*te(kt));
end

% combine magnitude weighted by normalised weights
combinedMagnitude = sum(magn .* (wTE ./ sum(wTE,4)),4);

end