function C = polyCondCheck(tdims, X, A, B, talk)
% Utility function returns the condition # estimate when A is preconditioned by a neumann series 
% composed of powers of B with coefficients given in X.
%>> tdims  Size of the coefficient matrix to use (if A and B not given)        int [3x1]
%>> X      Vector of series coefficients                                       double [Nx1]
%>> A      Operator to precondition (presumably a Laplacian)
%>> B      Off-diagonal parts after normalizing A
%>> talk   If existant, prints results after; Otherwise silent.

% If not given, create A and B operators
if nargin < 4
    A = createLaplacianMatrix(tdims);
    B = createBMatrix(tdims);
end

Bpoly = B;

% Identity operator
psi = X(1)*speye(prod(tdims));

% Compute polynomial sequence
for q = 1:(numel(X)-1);
    psi = psi + X(q+1)*Bpoly;

    Bpoly = Bpoly * B;
end

C = condest(psi*A);
if (nargin == 5) || (nargin == 3)
    fprintf('Cond of A: %g; Cond of P(x)A: %g\n', condest(A), C);
end

end
