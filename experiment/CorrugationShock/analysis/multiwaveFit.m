function [k a r] = multiwaveFit(q, dx, kguess, aguess)

vec = mkVector(kguess, aguess);

%r = waveIntegral(q, dx, vec);

for testloop = 1:30
%    H = computeHessian(q, dx, vec);
    F = computeGradient(q, dx, vec);

 %   vec = vec - (H^-1)*F;
    vec = vec - .25*F/norm(F);

%    disp([waveIntegral(q, dx, vec)]);
end

[k a] = unVector(vec);

r = waveIntegral(q, dx, vec);

end

% Computes the derivative of the functional with respect to each parameter
function F = computeGradient(q, dx, v)

epsilon = 1e-6;

F = zeros(size(v));

for j = 1:numel(v)
    vp = v;
    vp(j) = vp(j) + epsilon;
    f(1) = waveIntegral(q, dx, vp);
    
    vp = v;
    vp(j) = vp(j) - epsilon;
    f(2) = waveIntegral(q, dx, vp);

    F(j) = (f(1) - f(2))/(2*epsilon);
end

end

% Directly computes the matrix of second partials
function H = computeHessian(q, dx, v)

epsilon = 1e-6;

H = zeros(numel(v));

for j = 1:numel(v);
for k = j:numel(v);
    vp = v; vp(j) = vp(j) + epsilon; vp(k) = vp(k) + epsilon;
    f(1) = waveIntegral(q, dx, vp);

    vp = v; vp(j) = vp(j) - epsilon; vp(k) = vp(k) + epsilon;
    f(2) = waveIntegral(q, dx, vp);

    vp = v; vp(j) = vp(j) + epsilon; vp(k) = vp(k) - epsilon;
    f(3) = waveIntegral(q, dx, vp);

    vp = v; vp(j) = vp(j) - epsilon; vp(k) = vp(k) - epsilon;
    f(4) = waveIntegral(q, dx, vp);
    
    % Use the fact that the Hessian is symmetric for nice functions to save some time
    H(j,k) = (f(1) - f(2)) - (f(3) - f(4));
    H(k,j) = (f(1) - f(2)) - (f(3) - f(4));
end; end

H = H / (4*epsilon^2);

end

function v = mkVector(k, a);

if numel(k) ~= numel(a); error('Number of wavevectors and amplitudes not equal.'); end

v = [];

for j = 1:numel(k)
   v(end+1) = real(k(j));
   v(end+1) = imag(k(j));
end

for j = 1:numel(a);
   v(end+1) = real(a(j));
   v(end+1) = imag(a(j));
end

v = v(:);

end

function [k a] = unVector(v)

% The number of waves present = 1/4 the number of real parameters
N = numel(v) / 4;

x=1:N;

    k(x) = v(2*x-1) + 1i*v(2*x);
    a(x) = v(2*N + 2*x-1) + 1i*v(2*N + 2*x);

end

% Numerically evaluates our functional for q, given 
function ints = waveIntegral(q, dx, vec)
% q: The quantity to integrate
% x: The X values of the points in q
% t: The time values of the points in q
% a: An nx1 vector of wave amplitudes
% k: an nx1 vector of wavevectors
% w: a scalar, the frequency
% n: The number of waves to identify

[k a] = unVector(vec);
n = numel(k);

xvals = dx*(1:numel(q));
pred = zeros(size(q));

% Calculate the predicted perturbation from the given parameters
% All the variational quantities need this
for i = 1:n
	pred = pred + a(i)*exp(1i*(k(i)*xvals));
end

ints = sum(sum( (q - pred).*conj(q - pred) ));

end
