function [damprate corr wavevectorsY wavevectorsZ] = analyzeDampRates(dataframe, xran, yran, zran, dx)

sa = dataframe.mass(xran, :,:);

sr = zeros(size(sa));
for u = 1:size(sr, 1)
   sr(u,:,:) = fft2(squeeze(sa(u,:,:))); 
end

damprate = zeros(yran, zran);
corr = zeros(yran, zran);

% The x position of the cells
W = dx*(1:size(sr,1))';

for u = 1:yran; for v = 1:zran
    % Calculate d(log amp)/dx mode by mode -> exponential coefficents for exp. fit
    [f s]= polyfit(W, log(abs(sr(:,u,v))), 1);
    damprate(u,v) = f(1);
    corr(u,v) = s.normr;

end; end

wavevectorsY = (0:(yran-1))'*2*pi/((size(dataframe.mass,2)*dataframe.dGrid{2}));
wavevectorsZ = (0:(zran-1))'*2*pi/((size(dataframe.mass,2)*dataframe.dGrid{3}));

end

