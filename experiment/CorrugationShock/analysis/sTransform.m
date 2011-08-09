function [const res] = sTransform(X, u, dx, nwaves)
% X: the set of positions at which we know u.
% U: the vector of values to which we wish to fit a multiexponential curve
% dx: the cell spacing
% nwaves: How many terms in the multiexponential.

p0 = 10000;

for a = 2:numel(u)
	if u(a) < .5*p0; p0 = 1/X(a); break; end

	if a == numel(u); p0 = 1/X(round(end/2)); end % this case would suck
end

D = [];
Dpr = [];

for a = 0:(2*nwaves - 1)
    v = ((-X).^a) .* exp(-X*p0) .* u;
    D(a+1) = curveIntegral(v, X, dx) / factorial(a);
%   D(a+1) = sum(v)*dx / factorial(a);
end

%Dan = 1./[8 -8^2 8^3 -8^4 8^5 -8^6] + 1./[9 -9^2 9^3 -9^4 9^5 -9^6] + 1./[10 -10^2 10^3 -10^4 10^5 -10^6];

if nwaves == 1
	a0 = D(1);
	b1 = -D(2) / D(1);

	const = p0 - 1/b1;
	res = a0;
end

if nwaves == 2
	d0 = D(1); d1 = D(2); d2 = D(3); d3 = D(4);

	b1 = -(d1*d2 - d0*d3)/(d1*d1 - d2*d0);
	b2 = -(d1*d3 - d2*d2)/(d1*d1 - d2*d0);

	a0 = d0;
	a1 = b1*d0 + d1;

	q1 = (-b1 + sqrt(b1*b1-4*b2))/(2*b2);
	q2 = (-b1 - sqrt(b1*b1-4*b2))/(2*b2);

	const(1) = q1 + p0; const(2) = q2 + p0;
	res(1) = a0 + a1*const(1); res(2) = a0 + a1*const(2);
end

if nwaves == 3;
    d0 = D(1); d1 = D(2); d2 = D(3); d3 = D(4); d4 = D(5); d5 = D(6);

    B = -[d2 d1 d0; d3 d2 d1; d4 d3 d2]^-1 * [d3; d4; d5];
    A(1) = d0;
    A(2) = d1 + B(1)*d0;
    A(3) = d2 + B(1)*d1 + B(2)*d0;
    
    q = roots([B(3) B(2) B(1) 1]);
    
    const = q + p0;
    res = A(1) + A(2).*const + A(3) .* const.^2;
end

if nwaves == 4
    d0 = D(1); d1 = D(2); d2 = D(3); d3 = D(4); d4 = D(5); d5 = D(6); d6 = D(7); d7 = D(8);

    B = -[d3 d2 d1 d0; d4 d3 d2 d1; d5 d4 d3 d2; d6 d5 d4 d3]^-1 * [d4; d5; d6; d7];
    A(1) = d0;
    A(2) = d1 + B(1)*d0;
    A(3) = d2 + B(1)*d1 + B(2)*d0;
    A(4) = d3 + B(1)*d2 + B(2)*d1 + B(3);

    q = roots([B(4) B(3) B(2) B(1) 1]);

    const = q + p0;
    res = A(1) + A(2)*const + A(3)*const.^2 + A(4)*const.^3;
end

res = res*2*pi;

end
