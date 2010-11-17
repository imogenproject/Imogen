% Linear fit of time series
function [beta r] = linearValueFit(time, depvar,N )
beta = (sum(time .* depvar, 1) - (1/N)*sum(time,1).*sum(depvar,1)) ./ ...
       ( sum(time.^2,1) - (1/N)*sum(time,1).^2);

r = (N*sum(time.*depvar,1) - sum(time,1).*sum(depvar,1)) ./ (sqrt(N*sum(time.^2,1) - sum(time,1).^2).*sqrt(N*sum(depvar.^2,1) - sum(depvar,1).^2) );
end

