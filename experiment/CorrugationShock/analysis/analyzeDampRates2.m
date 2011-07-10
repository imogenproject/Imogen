function [kx omega kxRes omegaRes] = analyzeDampRates2(dq, x, t, thresh)
% > dq: Perturbed quantity
% > xvals: x position values for dq(ymode#, z mode#, x value)
yran = size(dq,1);
zran = size(dq,2);

damprate = zeros(size(dq,1), size(dq,2));
corr = zeros(size(damprate));

for u = 1:yran; for v = 1:zran
    % Calculate d(log amp)/dx mode by mode -> exponential coefficents for exp. fit
    Z = squeeze(dq(u,v,:,:));

    fitx = zeros(numel(t),2);
    sx = zeros(numel(t),1);
    usable = zeros(numel(t),1);

    for tp = 1:numel(t)
        Zcurrent = Z(:,tp);
        valids = (abs(Zcurrent) > thresh); % Select points whose amplitude is acceptable
        if numel(find(valids)) < 5; continue; end % Stop if we run out
        usable(tp) = 1;

        [alpha beta] = polyfit(x(valids), log(abs(Zcurrent(valids))), 1); % Make linear fit

        fitx(tp,:) = alpha(:);
        sx(tp) = beta.normr;
    end

    usable = (usable == 1);

    offsets = fitx(:,2); % Predicted dq at x=0 (at shock front)
    
    timefit = t(usable);
    offsets = offsets(usable);

    [wfit ws]= polyfit(timefit, offsets', 1);
    % Growth rate of offsets gives omega_im

    omega(u,v) = wfit(1);
    omegaRes(u,v) = ws.normr;

    % Predict kx by weighting by inverse residual.
    normalizer = [];
    if numel(find(usable == 1)) > 0;
        kxout = fitx(:,1);
        kxout = sum(kxout(usable)./abs(sx(usable)));

        normalizer = sum(1./abs(sx(usable)));

        kx(u,v) = kxout / normalizer;
        kxRes(u,v) = normalizer / numel(sx);
    else
        kx(u,v) = NaN;
        kxRes(u,v) = NaN;
    end
    fprintf('*');
end; end
fprintf('\n');

end

