% THIS IS Mathwork's bicgstab file, modified to work on the GPU.


function [x,flag,relres,iter,resvec] = bicgstab_GPU(A,b,tol,maxit,x0,varargin)
%BICGSTAB   BiConjugate Gradients Stabilized Method.
%   X = BICGSTAB(A,B) attempts to solve the system of linear equations
%   A*X=B for X. The N-by-N coefficient matrix A must be square and the
%   right hand side column vector B must have length N.
%
%   X = BICGSTAB(AFUN,B) accepts a function handle AFUN instead of the
%   matrix A. AFUN(X) accepts a vector input X and returns the
%   matrix-vector product A*X. In all of the following syntaxes, you can
%   replace A by AFUN.
%
%   X = BICGSTAB(A,B,TOL) specifies the tolerance of the method. If TOL is
%   [] then BICGSTAB uses the default, 1e-6.
%
%   X = BICGSTAB(A,B,TOL,MAXIT) specifies the maximum number of iterations.
%   If MAXIT is [] then BICGSTAB uses the default, min(N,20).
%
%   X = BICGSTAB(A,B,TOL,MAXIT,M) and X = BICGSTAB(A,B,TOL,MAXIT,M1,M2) use
%   preconditioner M or M=M1*M2 and effectively solve the system
%   A*inv(M)*X = B for X. If M is [] then a preconditioner is not
%   applied. M may be a function handle returning M\X.
%
%   X = BICGSTAB(A,B,TOL,MAXIT,M1,M2,X0) specifies the initial guess.  If
%   X0 is [] then BICGSTAB uses the default, an all zero vector.
%
%   [X,FLAG] = BICGSTAB(A,B,...) also returns a convergence FLAG:
%    0 BICGSTAB converged to the desired tolerance TOL within MAXIT iterations.
%    1 BICGSTAB iterated MAXIT times but did not converge.
%    2 preconditioner M was ill-conditioned.
%    3 BICGSTAB stagnated (two consecutive iterates were the same).
%    4 one of the scalar quantities calculated during BICGSTAB became
%      too small or too large to continue computing.
%
%   [X,FLAG,RELRES] = BICGSTAB(A,B,...) also returns the relative residual
%   NORM(B-A*X)/NORM(B). If FLAG is 0, then RELRES <= TOL.
%
%   [X,FLAG,RELRES,ITER] = BICGSTAB(A,B,...) also returns the iteration
%   number at which X was computed: 0 <= ITER <= MAXIT. ITER may be an
%   integer + 0.5, indicating convergence half way through an iteration.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = BICGSTAB(A,B,...) also returns a vector
%   of the residual norms at each half iteration, including NORM(B-A*X0).
%
%   Example:
%      n = 21; A = gallery('wilk',n);  b = sum(A,2);
%      tol = 1e-12;  maxit = 15; M = diag([10:-1:1 1 1:10]);
%      x = bicgstab(A,b,tol,maxit,M);
%   Or, use this matrix-vector product function
%      %-----------------------------------------------------------------%
%      function y = afun(x,n)
%      y = [0; x(1:n-1)] + [((n-1)/2:-1:0)'; (1:(n-1)/2)'].*x+[x(2:n); 0];
%      %-----------------------------------------------------------------%
%   and this preconditioner backsolve function
%      %------------------------------------------%
%      function y = mfun(r,n)
%      y = r ./ [((n-1)/2:-1:1)'; 1; (1:(n-1)/2)'];
%      %------------------------------------------%
%   as inputs to BICGSTAB:
%      x1 = bicgstab(@(x)afun(x,n),b,tol,maxit,@(x)mfun(x,n));
%
%   Class support for inputs A,B,M1,M2,X0 and the output of AFUN:
%      float: double
%
%   See also BICG, BICGSTABL, CGS, GMRES, LSQR, MINRES, PCG, QMR, SYMMLQ,
%   TFQMR, ILU, FUNCTION_HANDLE.

%   Copyright 1984-2008 The MathWorks, Inc.
%   $Revision: 1.19.4.4 $ $Date: 2008/10/31 06:20:38 $

% Check for an acceptable number of input arguments
if nargin < 2
    error('MATLAB:bicgstab:NotEnoughInputs', 'Not enough input arguments.');
end

    m = size(b,1);
    n = m;
    if ~isvector(b) || (size(b,2) ~= 1)
        error('MATLAB:bicgstab:RSHnotColumn',...
            'Right hand side must be a column vector.');
    end


% Assign default values to unspecified parameters
if nargin < 3 || isempty(tol)
    tol = 1e-6;
end
warned = 0;
if tol < eps
    warning('MATLAB:bicgstab:tooSmallTolerance', ...
        strcat('Input tol is smaller than eps and may not be achieved',...
        ' by BICGSTAB\n','         Try to use a bigger tolerance'));
    warned = 1;
    tol = eps;
elseif tol >= 1
    warning('MATLAB:bicgstab:tooBigTolerance', ...
        strcat('Input tol is bigger than 1 \n',...
        '         Try to use a smaller tolerance'));
    warned = 1;
    tol = 1-eps;
end
if nargin < 4 || isempty(maxit)
    maxit = min(n,20);
end

% Check for all zero right hand side vector => all zero solution
n2b = gpunorm(b);                     % Norm of rhs vector, b
if (n2b == 0)                      % if    rhs vector is all zeros
    x = zeros(n,1);                % then  solution is all zeros
    flag = 0;                      % a valid solution has been obtained
    relres = 0;                    % the relative residual is actually 0/0
    iter = 0;                      % no iterations need be performed
    resvec = 0;                    % resvec(1) = gpunorm(b-A*x) = gpunorm(0)
    if (nargout < 2)
        itermsg('bicgstab',tol,maxit,0,flag,iter,NaN);
    end
    return
end

%if ((nargin >= 5) && ~isempty(M1))
%    existM1 = 1;
%    [m1type,m1fun,m1fcnstr] = iterchk(M1);
%    if strcmp(m1type,'matrix')
%        if ~isequal(size(M1),[m,m])
%            error('MATLAB:bicgstab:WrongPrecondSize', ...
%                ['Preconditioner must be a square matrix' ...
%                ' of size %d to match the problem size.'],m);
%        end
%    end
%else
%    existM1 = 0;
%    m1type = 'matrix';
%end
%
%if ((nargin >= 6) && ~isempty(M2))
%    existM2 = 1;
%    [m2type,m2fun,m2fcnstr] = iterchk(M2);
%    if strcmp(m2type,'matrix')
%        if ~isequal(size(M2),[m,m])
%            error('MATLAB:bicgstab:WrongPrecondSize', ...
%                ['Preconditioner must be a square matrix' ...
%                ' of size %d to match the problem size.'],m);
%        end
%    end
%else
%    existM2 = 0;
%%    m2type = 'matrix';
%end
%
if ((nargin >= 5) && ~isempty(x0))
    if ~isequal(size(x0),[n,1])
        error('MATLAB:bicgstab:WrongInitGuessSize', ...
            ['Initial guess must be a column vector of' ...
            ' length %d to match the problem size.'],n);
    else
        x = x0;
    end
else
    x = GPUdouble(zeros(n,1));
end

if ((nargin > 5) && strcmp(atype,'matrix') && ...
        strcmp(m1type,'matrix') && strcmp(m2type,'matrix'))
    error('MATLAB:bicgstab:TooManyInputs', 'Too many input arguments.');
end


% Set up for the method
flag = 1;
xmin = x;                          % Iterate which has minimal residual so far
imin = 0;                          % Iteration at which xmin was computed
tolb = tol * n2b;                  % Relative tolerance
r = b - A(x);
%%%%r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
normr = gpunorm(r);                   % Norm of residual
normr_act = normr;

if (normr <= tolb)                 % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2b;
    iter = 0;
    resvec = normr;
    if (nargout < 2)
        itermsg('bicgstab',tol,maxit,0,flag,iter,relres);
    end
    return
end

rt = r;                            % Shadow residual
resvec = zeros(2*maxit+1,1);       % Preallocate vector for norm of residuals
resvec(1) = normr;                 % resvec(1) = gpunorm(b-A*x0)
normrmin = normr;                  % Norm of residual from xmin
rho = 1;
omega = 1;
stag = 0;                          % stagnation of the method
alpha = [];                        % overshadow any functions named alpha
moresteps = 0;
maxmsteps = min([floor(n/50),10,n-maxit]);
maxstagsteps = 3;
% loop over maxit iterations (unless convergence or failure)

for ii = 1 : maxit
    rho1 = rho;
    %rho = sum(rt .* r);
    rho = dotProd(rt, r);

%    rho = rt' * r;
    if (rho == 0.0) || isinf(rho)
        flag = 4;
        resvec = resvec(1:2*ii-1);
        break
    end
    if ii == 1
        p = r;
    else
        beta = (rho/rho1)*(alpha/omega);
        if (beta == 0) || ~isfinite(beta)
            flag = 4;
            break
        end
        p = r + beta * (p - omega * v);
    end

    v = A(p);
    rtv = dotProd(rt, v);
%    rtv = sum(rt .* v);
%    rtv = rt' * v;
    
    if (rtv == 0) || isinf(rtv)
        flag = 4;
        resvec = resvec(1:2*ii-1);
        break
    end
    alpha = rho / rtv;
    if isinf(alpha)
        flag = 4;
        resvec = resvec(1:2*ii-1);
        break
    end
    
    if abs(alpha)*gpunorm(p) < eps*gpunorm(x)
        stag = stag + 1;
    else
        stag = 0;
    end
    
    xhalf = x + alpha * p;        % form the "half" iterate
    s = r - alpha * v;             % residual associated with xhalf
    normr = gpunorm(s);
    normr_act = normr;
    resvec(2*ii) = normr;
    
    % check for convergence
    if (normr <= tolb || stag >= maxstagsteps || moresteps)
        s = b - A(xhalf);

        normr_act = gpunorm(s);
        resvec(2*ii) = normr_act;
        if normr_act <= tolb
            x = xhalf;
            flag = 0;
            iter = ii - 0.5;
            resvec = resvec(1:2*ii);
            break
        else
            if stag >= maxstagsteps && moresteps == 0
                stag = 0;
            end
            moresteps = moresteps + 1;
            if moresteps >= maxmsteps
                if ~warned
                    warning('MATLAB:bicgstab:tooSmallTolerance', ...
                        strcat('Input tol may be smaller than eps*cond(A)',...
                        ' and might not be achieved by BICGSTAB\n',...
                        '         Try to use a bigger tolerance'));
                end
                flag = 3;
                x = xhalf;
                resvec = resvec(1:2*ii);
                break;
            end
        end
    end
    
    if stag >= maxstagsteps
        flag = 3;
        resvec = resvec(1:2*ii);
        break
    end
    
    if normr_act < normrmin        % update minimal norm quantities
        normrmin = normr_act;
        xmin = xhalf;
        imin = ii - 0.5;
    end
    
    t = A(s);
 
    tt = dotProd(t, t);   
%    tt = sum(t .* t);

    if (tt == 0) || isinf(tt)
        flag = 4;
        resvec = resvec(1:2*ii);
        break
    end
    omega = dotProd(t, s) / tt;
%    omega = (sum(t .* s)) / tt;
    if isinf(omega)
        flag = 4;
        resvec = resvec(1:2*ii);
        break
    end
    
    if abs(omega)*gpunorm(s) < eps*gpunorm(xhalf)
        stag = stag + 1;
    else
        stag = 0;
    end
    
    x = xhalf + omega * s;        % x = (x + alpha * ph) + omega * sh
    r = s - omega * t;
    normr = gpunorm(r);
    normr_act = normr;
    resvec(2*ii+1) = normr;
    
    % check for convergence        
    if (normr <= tolb || stag >= maxstagsteps || moresteps)
        r = b - A(x);
%        r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
        normr_act = gpunorm(r);
        resvec(2*ii+1) = normr_act;
        if normr_act <= tolb
            flag = 0;
            iter = ii;
            resvec = resvec(1:2*ii+1);
            break
        else
            if stag >= maxstagsteps && moresteps == 0
                stag = 0;
            end
            moresteps = moresteps + 1;
            if moresteps >= maxmsteps
                if ~warned
                    warning('MATLAB:bicgstab:tooSmallTolerance', ...
                        strcat('Input tol may be smaller than eps*cond(A)',...
                        ' and might not be achieved by BICGSTAB\n',...
                        '         Try to use a bigger tolerance'));
                end
                flag = 3;
                resvec = resvec(1:2*ii+1);
                break;
            end
        end        
    end
    
    if normr_act < normrmin        % update minimal norm quantities
        normrmin = normr_act;
        xmin = x;
        imin = ii;
    end
    
    if stag >= maxstagsteps
        flag = 3;
        resvec = resvec(1:2*ii+1);
        break
    end
    
end                                % for ii = 1 : maxit

% returned solution is first with minimal residual
if flag == 0
    relres = normr_act / n2b;
else
    r = b - A(xmin);
%    r = b - iterapp('mtimes',afun,atype,afcnstr,xmin,varargin{:});
    if gpunorm(r) <= normr_act
        x = xmin;
        iter = imin;
        relres = gpunorm(r) / n2b;
    else
        iter = ii;
        relres = normr_act / n2b;
    end
end

% only display a message if the output flag is not used
if nargout < 2
    itermsg('bicgstab',tol,maxit,ii,flag,iter,relres);
end

end

function D = dotProd(a, b)
D = wrap_cublasDdot(a,b);
end

function n = gpunorm(x)
n = sqrt(wrap_cublasDdot(x,x));
end
