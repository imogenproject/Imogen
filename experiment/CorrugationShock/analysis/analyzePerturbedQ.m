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
    kx_im_residual = zeros(numel(x),1);
    fitw_re = zeros(numel(x),2);
    kx_re_residual = zeros(numel(x),1);
    
    usable = zeros(numel(t),1);

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

    % Use constant offsets for omega to predict amplitude at constant T,
    % thereby revealing f(x) and assume single exponential decay.
    
    % %%%%%%%%%%%%%% Fit kx decay
    offsets = fitw_im(:,2); % Predicted |dq| at t=0

    spacefit = x(usable);
    offsets = offsets(usable);

    [xfit ws]= polyfit(spacefit, offsets, 1);
    % Growth rate of offsets gives kx_im

    kxImag = xfit(1);
    kxImagRes = ws.normr;
    
    %%%%%%%%%% Fit kx oscillatory
    offsets = fitw_re(:,2); % Predicted phase(dq) at t=0
    
    spacefit = x(usable);
    offsets = offsets(usable);

    [xfit ws]= polyfit(spacefit, offsets, 1);
    % Growth rate of offsets gives kx_im

    kxRe = xfit(1);
    kxReRes = ws.normr;

    kx(u,v) = kxRe + 1i*kxImag;
    kxRes(u,v) = kxReRes + 1i*kxImagRes;
    
    % Predict w by weighting by inverse residual.
    normalizer = [];
    if numel(find(usable == 1)) > 0;
        % Fit decay rates
        wout = fitw_im(:,1);
        wout = sum(wout(usable)./abs(w_im_residual(usable))');

        normalizer = sum(1./abs(w_im_residual(usable)));

        omega(u,v) = 1i*wout / normalizer;
        omegaRes(u,v) = 1i* normalizer / numel(w_im_residual);

        % Fit oscillation rates
        wout = fitw_re(:,1);
        wout = sum(wout(usable)./abs(w_re_residual(usable))');

        normalizer = sum(1./abs(w_re_residual(usable)));

        omega(u,v) = omega(u,v) + wout / normalizer;
        omegaRes(u,v) = omegaRes(u,v) + normalizer / numel(w_im_residual);
    else
        omega(u,v) = NaN;
        omegaRes(u,v) = NaN;
    end
    
    fprintf('*');
end; end
fprintf('\n');

end

