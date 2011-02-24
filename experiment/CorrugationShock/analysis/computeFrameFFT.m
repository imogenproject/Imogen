function result = computeFrameFFT(frame, xset, yset, zset)

%dim = size(frame.mass);
%xhalf = round(dim(1)/2);

%fourierSet = frame.momY(xset, yset, zset) ./ frame.mass(xset, yset, zset);
fourierSet = frame.mass(xset, yset, zset);

result = zeros(size(fourierSet));

for n = 1:size(result,1)
	result(n,:,:) = fft2(squeeze(fourierSet(n,:,:)));
end

end
