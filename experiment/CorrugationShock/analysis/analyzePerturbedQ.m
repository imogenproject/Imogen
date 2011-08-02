function [kx omega kxOffset omegaOffset kxRes omegaRes] = analyzePeturbedQ(dq, x, t, linearFrames, preorpost)
% > dq: Perturbed quantity
% > xvals: x position values for dq(ymode#, z mode#, x value)
yran = size(dq,1);
zran = size(dq,2);

for u = 1:yran; for v = 1:zran

    % Omega fit
    if strcmp(preorpost,'post')
        [wimfit confidenceIm] = polyfit(t(linearFrames), mean(squeeze(log(abs(dq(u,v,2:15,linearFrames))))),1);
        [wrefit confidenceRe] = polyfit(t(linearFrames), mean(unwrap(squeeze(angle(dq(u,v,2:15,linearFrames))),1,2 )), 1);

        [kximfit confidenceKxIm] = polyfit(x, mean(squeeze(log(abs(dq(u,v,:,linearFrames((end-10):end))))),2),1);
        [kxrefit confidenceKxRe] = polyfit(x, mean(unwrap(squeeze(angle(dq(u,v,:,linearFrames((end-10):end)))),pi,1),2), 1);
    else
        [wimfit confidenceIm] = polyfit(t(linearFrames), mean(squeeze(log(abs(dq(u,v,(end-10):(end-1),linearFrames))))),1);
        [wrefit confidenceRe] = polyfit(t(linearFrames), mean(unwrap(squeeze(angle(dq(u,v,(end-10):(end-1),linearFrames))),1,2 )), 1);

        [kximfit confidenceKxIm] = polyfit(x, mean(squeeze(log(abs(dq(u,v,:,linearFrames((end-10):end))))),2),1);
        [kxrefit confidenceKxRe] = polyfit(x, mean(unwrap(squeeze(angle(dq(u,v,:,linearFrames((end-10):end)))),pi,1),2), 1);
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

