function result = computeFrameFFT(frame, xset, yset, zset)

%dim = size(frame.mass);
%xhalf = round(dim(1)/2);

fourierSet = frame.momY(xset, yset, zset) ./ frame.mass(xset, yset, zset);
result = fft2(squeeze(fourierSet));

end
