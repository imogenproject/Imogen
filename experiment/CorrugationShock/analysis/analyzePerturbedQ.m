function [kx omega kxOffset omegaOffset kxRes omegaRes] = analyzePeturbedQ(dq, x, t, linearFrames, preorpost)
% > dq: Perturbed quantity
% > xvals: x position values for dq(ymode#, z mode#, x value)
yran = size(dq,1);
zran = size(dq,2);

for u = 1:yran; for v = 1:zran

    % Omega fit
    if strcmp(preorpost,'post')
        
        [wimfit confidenceIm] = monovariateFit(t(linearFrames), mean(squeeze(log(abs(dq(u,v,2:15,linearFrames))))));
        [wrefit confidenceRe] = monovariateFit(t(linearFrames), mean(unwrap(squeeze(angle(dq(u,v,2:15,linearFrames))),1,2 )));

        [kximfit confidenceKxIm] = monovariateFit(x, mean(squeeze(log(abs(dq(u,v,:,linearFrames((end-10):end))))),2));
        [kxrefit confidenceKxRe] = monovariateFit(x, mean(unwrap(squeeze(angle(dq(u,v,:,linearFrames((end-10):end)))),pi,1),2));
    else
        [wimfit confidenceIm] = monovariateFit(t(linearFrames), mean(squeeze(log(abs(dq(u,v,(end-10):(end-1),linearFrames))))));
        [wrefit confidenceRe] = monovariateFit(t(linearFrames), mean(unwrap(squeeze(angle(dq(u,v,(end-10):(end-1),linearFrames))),1,2 )));

        [kximfit confidenceKxIm] = monovariateFit(x, mean(squeeze(log(abs(dq(u,v,:,linearFrames((end-10):end))))),2));
        [kxrefit confidenceKxRe] = monovariateFit(x, mean(unwrap(squeeze(angle(dq(u,v,:,linearFrames((end-10):end)))),pi,1),2));
    end
    omega(u,v) = wrefit(1) + 1i*wimfit(1);
    omegaRes(u,v) = confidenceRe.normr + 1i*confidenceIm.normr;
    omegaOffset(u,v) = wrefit(2) + 1i*wimfit(2);

    kx(u,v) = kxrefit(1) + 1i*kximfit(1);
    kxRes(u,v) = confidenceKxRe.normr + 1i*confidenceKxIm.normr;
    kxOffset(u,v) = kxrefit(2) + 1i*kximfit(2);
 
    fprintf('*');
end; end
fprintf('\n');

end

function [fit, residual] = monovariateFit(x, y)
%[fit, residual] = polyfit(x, y, 1);

x = x(:); y = y(:); %make it nx1

%deriv2 = circshift(y,-1) + circshift(y,1) - 2*y; % Calculate second derivative

%deriv2(1) = deriv2(2);
%deriv2(end) = deriv2(end-1);

%dx = x - circshift(x,1); dx(end) = dx(end-1);

%deriv2 = deriv2 ./ dx;

%deriv2 = deriv2 / max(deriv2);

%weight = (1 + 2*deriv2).^-1; % FIXME: parameterize this
%weight = weight / sum(weight);

weight = ones(size(x));

N = numel(x);

ans = ([N sum(x.*weight); sum(x.*weight) sum(x.^2 .* weight)]^-1) * [sum(y.*weight); sum(x.*y.*weight)];

fit = [ans(2) ans(1)];

residual.normr = sqrt(sum(y - (fit(2) + fit(1)*x)));

end
