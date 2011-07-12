function [kx omega kxRes omegaRes] = analyzePeturbedQ(dq, x, t, thresh)
% > dq: Perturbed quantity
% > xvals: x position values for dq(ymode#, z mode#, x value)
yran = size(dq,1);
zran = size(dq,2);

damprate = zeros(size(dq,1), size(dq,2));
corr = zeros(size(damprate));

for u = 1:yran; for v = 1:zran
    % Calculate d(log amp)/dx mode by mode -> exponential coefficents for exp. fit
    Z = squeeze(dq(u,v,:,:));

    fitkx_im = zeros(numel(t),2);
    kx_im_residual = zeros(numel(t),1);
    fitkx_re = zeros(numel(t),2);
    kx_re_residual = zeros(numel(t),1);
    
    usable = zeros(numel(t),1);

    for tp = 1:numel(t)
        Zcurrent = Z(:,tp);
        valids = (abs(Zcurrent) > thresh); % Select points whose amplitude is acceptable
        if numel(find(valids)) < 5; continue; end % Stop if we run out
        usable(tp) = 1;

        [alpha beta] = polyfit(x(valids), log(abs(Zcurrent(valids))), 1); % Make linear fit
        fitkx_im(tp,:) = alpha(:);
        kx_im_residual(tp) = beta.normr;

        [alpha beta] = polyfit(x(valids), unwrap(angle(Zcurrent(valids))),1); % linear fit phase(x)
        fitkx_re(tp,:) = alpha(:);
        kx_re_residual(tp) = beta.normr;
    end

    usable = (usable == 1);

    offsets = fitkx_im(:,2); % Predicted |dq| at x=0 (at shock front)
    
    timefit = t(usable);
    offsets = offsets(usable);

    [wfit ws]= polyfit(timefit, offsets', 1);
    % Growth rate of offsets gives omega_im

    omegaImag = wfit(1);
    omegaImagRes = ws.normr;

    % Predict kx by weighting by inverse residual.
    normalizer = [];
    if numel(find(usable == 1)) > 0;
        % Fit decay rates
        kxout = fitkx_im(:,1);
        kxout = sum(kxout(usable)./abs(kx_im_residual(usable)));

        normalizer = sum(1./abs(kx_im_residual(usable)));

        kx(u,v) = 1i*kxout / normalizer;
        kxRes(u,v) = 1i* normalizer / numel(kx_im_residual);

        % Fit oscillation rates
        kxout = fitkx_re(:,1);
        kxout = sum(kxout(usable)./abs(kx_re_residual(usable)));

        normalizer = sum(1./abs(kx_re_residual(usable)));

        kx(u,v) = kx(u,v) + kxout / normalizer;
        kxRes(u,v) = kxRes(u,v) + normalizer / numel(kx_im_residual);
    else
        kx(u,v) = NaN;
        kxRes(u,v) = NaN;
    end
    
    offsets = unwrap(fitkx_re(usable,2)); % predicted phase at |x|=0

    [wfit ws]= polyfit(timefit, offsets', 1);
    % phase rate of offsets gives us omega_re

    omega(u,v) = wfit(1) + 1i*omegaImag;
    omegaRes(u,v) = ws.normr + 1i*omegaImagRes;

    fprintf('*');
end; end
fprintf('\n');

end

