function Ax = polyPreconditionL2_GPU(x, dims, noIniReshape)
% Polynomial preconditioner for Laplacian operators

% Polynomials with 1-6 terms are effective. Polynomials with more do not help sufficiently.
numPolyTerms = 7;
polycoeffs = ...
[ 1      0      0      0      0      0      0      0      0      0      0      0; ...
  .86666 1.1666 0      0      0      0      0      0      0      0      0      0; ...
  1.1430 1.6130 1.0930 0      0      0      0      0      0      0      0      0; ...
  0.6148 0.2997 0.3457 4.5327 0      0      0      0      0      0      0      0; ...
...%  0.9250 1.3250 2.2750 1.625  0      0      0      0      0      0      0      0; ...
  0.6637 0.5379 0      1.4652 4.6500 0      0      0      0      0      0      0 ; ...
...%  0.6989 0.0694 0.4155 2.4365 2.1069 0      0      0      0      0      0      0; ...
  0.7035 0.6501 0.0506 0      3.0529 4.8318 0 0 0 0 0 0; ...
...%  0.7831 0.6615 0.0940 1.0152 3.0936 4.2213 0 0 0 0 0 0; ...
...%  1.01   0.5046 0.7589 0.9939 2.9454 2.9431 0      0      0      0      0      0; ...
  0.6134 0.5883 0.1932 0.1571 1.0606 2.8089 4.1310 0 0 0 0 0; ...
...%  1      0.3976 0.6940 1.2421 1.5816 3.0061 2.2642 0      0      0      0      0; ...
  0.7507 0.5168 0.7383 0.3834 0.8655 0.9014 1.7078 2.1066 0      0      0      0; ...
  0.8382 0.4303 0.6870 0.7487 0.6205 1.1311 1.8124 2.3457 1.4084 0      0      0; ...
  0.9777 0.6226 0.9598 0.7665 0.9870 0.8902 1.2504 1.4019 1.2811 1.5306 0      0; ...
  0.8985 0.5576 0.6832 0.6565 1.2490 0.7570 0.7062 0.9389 0.8090 2.5012 1.5101 0; ...
  0.6152 0.3859 0.4257 0.4166 0.2121 0.1266 1.2330 0.8020 0.5978 1.8693 2.0348 1.2110];

if ~noIniReshape
    x   = reshape(x, dims);
end

Npoints = prod(dims);

% Let Ax denote the entire preconditioning operator acting on x
% Let Bn_x denote the nth power of the B operator acting on x
%--- Store first term [identity operator * x] in Ax; Prealloc q ---%
Ax  = GPUdouble(); setReal(Ax); setSize(Ax, dims); GPUallocVector(Ax);
GPUtimes(x, polycoeffs(numPolyTerms,1), Ax);
%Ax = polycoeffs(numPolyTerms,1) * x;

Bn_x = GPUdouble(); setReal(Bn_x); setSize(Bn_x, dims); GPUallocVector(Bn_x);
cublasScopy(2*Npoints, getPtr(x), 1, getPtr(Bn_x), 1);
%Bn_x = 1.0 * x;

accum = GPUdouble(); setReal(accum); setSize(accum, dims); GPUallocVector(accum);
accum(:)=0;

%--- Iterate up to requested polynomial order ---%
for ord = 2:numPolyTerms;
    %--- Perform term summation, alternating source & overwritten buffer ---%
    if mod(ord, 2) == 0;
        accumulateBterm(Bn_x, accum, polycoeffs(numPolyTerms, ord), Ax);
    else
        accumulateBterm(accum, Bn_x, polycoeffs(numPolyTerms, ord), Ax);
    end

end

clear Bn_x;
Ax = reshape(Ax, [Npoints 1]);

end
