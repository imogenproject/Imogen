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

    fitw_im = zeros(numel(x),2);
    w_im_residual = zeros(numel(x),1);
    fitw_re = zeros(numel(x),2);
    w_re_residual = zeros(numel(x),1);
    
    usable = zeros(numel(x),1);

    % Fit all curves of constant x to get many omega values
    for tp = 1:numel(x)
        Zcurrent = Z(tp,:);
        valids = (abs(Zcurrent) > thresh); % Select points whose amplitude is acceptable
        if numel(find(valids)) < 5; continue; end % Stop if we run out
        usable(tp) = 1;

        [alpha beta] = polyfit(t(valids), log(abs(Zcurrent(valids))), 1); % Make linear fit
        fitw_im(tp,:) = alpha(:);
        w_im_residual(tp) = beta.normr;

        [alpha beta] = polyfit(t(valids), unwrap(angle(Zcurrent(valids))),1); % linear fit phase(x)
        fitw_re(tp,:) = alpha(:);
        w_re_residual(tp) = beta.normr;
    end

    usable = (usable == 1);
    
    % Predict w by weighting by inverse residual.
    normalizer = [];
    if numel(find(usable == 1)) > 0;
        % Fit decay rates
        wout = fitw_im(:,1);
        wout = sum(wout(usable)./abs(w_im_residual(usable)));

        normalizer = sum(1./abs(w_im_residual(usable)));

        omega(u,v) = 1i*wout / normalizer;
        omegaRes(u,v) = 1i* normalizer / numel(w_im_residual);

        % Fit oscillation rates
        wout = fitw_re(:,1);
        wout = sum(wout(usable)./abs(w_re_residual(usable)));

        normalizer = sum(1./abs(w_re_residual(usable)));

        omega(u,v) = omega(u,v) + wout / normalizer;
        omegaRes(u,v) = omegaRes(u,v) + normalizer / numel(w_im_residual);
    else
        omega(u,v) = NaN;
        omegaRes(u,v) = NaN;
    end
    
    
    %%%%%%%%%%%%%%% fit kx
    fitkx_im = zeros(numel(t),2);
    kx_im_residual = zeros(numel(t),1);
    fitkx_re = zeros(numel(t),2);
    kx_re_residual = zeros(numel(t),1);
    
    usable = zeros(numel(t),1);

    % Fit all curves of constant t to get many kx values
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

        % Predict kx by weighting with inverse residual.
    normalizer = [];
    usable = (usable == 1);
    
    if numel(find(usable == 1)) > 0;
        % Fit decay rates
        kxout = fitkx_im(:,1);
        kxout = sum(kxout(usable)./abs(kx_im_residual(usable)));

        normalizer = sum(1./abs(kx_im_residual(usable)));

        kx(u,v) = 1i*wout / normalizer;
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
    
    
    fprintf('*');
end; end
fprintf('\n');

end

