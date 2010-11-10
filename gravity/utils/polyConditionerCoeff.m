function [X hist] = polyConditionerCoeff(tdims, imax, X, resumehist)
% Utility function attempts to search coefficient parameter space of polynomial preconditioner's
% Neumann series to identify the best coefficients. Employs gradient descent N-dimensional search.
%>>  tdims Size of the coefficient matrix to generate for the search           int [3x1]
%>>  imax  Maximum number of steps to take                                     int
%>>  x0    Optional initial position parameter                                 double [Nx1]
%<<  X     Final iterate with best coefficients                                double [Nx1]
%<<  hist  History of condition numbers with accompanying coeffcients     double [N+1 by imax+1]

% Step size
h = .01;
maxEscapeTries = 3;
nEscapeTries = 0;

normFunction = @polyCondCheck;

%--- Create needed linear operator matrices ---%
A = createCSOMatrix(tdims, [-24 2 1 0]);
B = createCSOMatrix(tdims, [0 2 1 0]/24);
ID = speye(prod(tdims));

keepsearching = 1;
niter = 0;

fprintf('Cond of P(Xini)A: %g\n', normFunction(tdims, X, A, B));

stepcoeff = .05 * ones(size(X));
stepMAX = size(X,1)/20;

if nargin == 4
    hist = resumehist;
    niter = size(hist,1) + 1;
else
    niter = 1;
end

%--- Initiate search ---%
%        Sample the condition at X. In parallel, compute forward derivatives for every dimension.
%        Attempt to descend gradient towards minimal condition number.
%        Print some progress info so people don't give up in utter hopeless despair.

% Store initial guess in return vector
hist(niter,1) = normFunction(tdims, X, A, B);
for z = 1:numel(X); hist(niter,1+z) = X(z); end

nmax = numel(X);
delta = zeros([1 nmax]);


while keepsearching
    niter = niter + 1;

    f0 = hist(niter-1,1);
    Xp = repmat(X, [nmax 1]);
    for z = 1:nmax; Xp(z,z) = Xp(z,z) + h; end

    delta2 = delta;
    parfor z = 1:nmax;
        fp = normFunction(tdims, Xp(z,:), A, B);
        delta(z) = -.1 * sign(fp-f0) * stepcoeff(z) / (h*max(abs(fp-f0), 10));
    end
    X = X + delta;
    X(X < 0) = 0;
    if X(1) < .5; X(1)=.5; end

    % If we're oscillating, shrink step rapidly. Otherwise grow towards a maximum
    for z = 1:nmax; if (sign(delta2(z)) + sign(delta(z))) == 0; stepcoeff(z) = stepcoeff(z) * .5; else; stepcoeff(z) = stepcoeff(z) * 1.05 / (1 + .05*(stepcoeff(z)/stepMAX)^2); end; end

    hist(niter,1) = normFunction(tdims, X, A, B);
    for z = 1:numel(X); hist(niter,1+z) = X(z); end

    % If the iterates appear to stall, attempt to escape local minimum
    % Don't do this if the norm rises since that may be because of jittering last iteration
    % If we are making progress
    if (hist(niter,1) - hist(niter-1,1)) < 0
        % But it's really slow
        if abs(hist(niter,1) - hist(niter-1,1)) < .005
            % Compare this to the best guess yet
            condBest = min(hist(:,1));

            % If this damn near the best guess yet, we only try so many times before accepting it and going into asymptotic approach
            if (abs(hist(niter,1) - condBest)) < .01
                nEscapeTries = nEscapeTries + 1;

                if nEscapeTries <= maxEscapeTries
                    fprintf('\nStalled at cond=%g near best of %g %i/%i; Jittering...\n', hist(niter,1), condBest, nEscapeTries, maxEscapeTries); 
                    X = doTheJitterbug(X, nEscapeTries);
                    stepcoeff = getDefaultStep(X);
               end

            else
                fprintf('\nStalled at cond=%g; Jittering...\n', hist(niter,1));
                X = doTheJitterbug(X, 0);
                stepcoeff = getDefaultStep(X);
            end

        end

        
    end
  
    if mod(niter, 10) == 0; fprintf('X'); else fprintf('*'); end
    if niter >= imax; keepsearching = 0; end
end

% Identify best iterate and store as return Z
[z IDX] = sort(hist(:,1));
z = IDX(1);

X = hist(z,2:end);

fprintf('\nCond of P(Xbest)A: %g\n', normFunction(tdims, X, A, B));

end

function U = doTheJitterbug(X, trynum)
U = X + sign(rand(size(X))-.5).*rand(size(X))*(.25+.1*trynum);
end

function st = getDefaultStep(n)
st = .05 * ones(size(n));
end
