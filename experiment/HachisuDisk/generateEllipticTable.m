function struct =  generateEllipticTable(delta)
% This function computes a table of complete elliptic integrals.
% The integral of 1/sqrt(a - b cos theta) can be transformed into
% 1/sqrt(a) * integral of 1/sqrt(1 - (b/a) cos theta ) for which an interpolation table requires
% only a single parameter. For cylindrical gravitational potential integrals this falls in the
% range for which the integral is real, [0 1);
rinv = @(angle, param) 1 ./ sqrt(1 - param*cos(angle));

BoverA = -delta:delta:(1+delta);
Fi = zeros(size(BoverA));

for alpha = 2:(size(BoverA,2)-2)
    Fi(alpha) = 2*quadl(@(theta) rinv(theta, BoverA(alpha)), 0, pi, 1e-9);
end

% Avoids problems with NaN and divergence when we get close to 1
Fi(end-1) = Fi(end-2);
Fi(end) = Fi(end-1);
% Makes for easy derivatives
Fi(1) = Fi(2);

struct.BoverA = BoverA;
struct.F = Fi;
struct.df_dq = (circshift(Fi, [0 -1]) - circshift(Fi, [0 1])) / (2*delta);

struct.df_dq(end-2) = (struct.F(end-2) - struct.F(end-3))/delta;
struct.df_dq(end-1:end) = struct.df_dq(end-2);

end
