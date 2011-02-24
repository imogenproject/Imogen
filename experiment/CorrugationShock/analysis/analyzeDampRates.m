function [damprate corr mdr] = findDampRates(dataframe, xran, yran, zran, fbar, dx)

sr = log(abs(dataframe.mass(xran, yran, zran) - fbar));

damprate = zeros(size(sr,2), size(sr,3) );
corr = zeros(size(sr,2), size(sr,3) );

% The x position of the cells
W = dx*(1:size(sr,1))';

for u = 1:size(sr,2); for v = 1:size(sr,3)
    % Calculate d(log amp)/dx -> exponential for exp. fit
    [f s]= polyfit(W, sr(:,u,v), 1);
    damprate(u,v) = f(1);
    corr(u,v) = s.normr;

end; end

% Mean Damp Rate - provide an inverse-residual-weighted mean estimate of the damp rate
mdr = sum(sum(damprate ./ corr)) / sum(sum(1./corr));

end

