function result = computeFrameFFT(frame)

dim = size(frame.mass);

xhalf = round(dim(1)/2);

rhohalf = frame.magZ(xhalf + 5,:,:);
result = fft2(squeeze(rhohalf));

end
