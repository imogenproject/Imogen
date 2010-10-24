% Function rudely prevents velocities from exceeding a given max.
% Maybe this will kill those sodding point-explosions in the disk simulations!

function [r px py e] = speedLimiter(rho, momx, momy, ener, vmax)

grid = size(momx);

vsqr = (momx.^2 + momy.^2) ./ (rho.^2);

toofast = (simpleBlur((vsqr > vmax^2), 1, 1) > 0);

I = find(toofast); % Grab the locations of all the too-fast points

if numel(I) == 0
r = rho; px = momx; py = momy; e = ener;
return;
end

fprintf('Speed limiter: Modified %i cells\n', numel(I));

f = simpleBlur(rho, 2, .5); r = rho;   r(I) = f(I);
f = simpleBlur(momx,2, .5); px = momx; px(I) = f(I);
f = simpleBlur(momy,2, .5); py = momy; py(I) = f(I);
f = simpleBlur(ener,2, .5); e = ener;  e(I) = f(I);

end
