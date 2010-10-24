function rOuter = kojimaFindOuterRadius(q, radiusRatio)
% This computes the outer radius of a Kojima disk given two parameters. Curiously radius
% does not depend on the polytropic index
%
%>> q           Angular velocity omega scales as R^(2-2*q); q = 1.5 -> Keplerian motion  double
%>> radiusRatio Radius of inner edge to radius of density max                            double
%<< rOuter      Outer radius of the disk at the equatorial plane                         double

xq = 2-2*q;
guess = radiusRatio / (2.0 * radiusRatio - 1.0);

%--- Newton's method; Given the above guess generally converges in < 4 iterations ---%
for loopVar = 1:20;
    f = (radiusRatio^xq-guess^xq)-xq*(1/guess-1/radiusRatio);
    df = xq * (guess^-2 - guess^(xq-1));

    delta = -f / df;

    if abs(delta/guess) <= 1e-8; rOuter = guess; return; end;

    guess = guess + delta;
end

rOuter = guess;

end
