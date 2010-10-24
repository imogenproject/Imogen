function phi = bicgstabPotentialSolver_GPU(run, mass, phi0)
% This routine solves for gravitational potential using a linear solver on the interior and 
% finding potential at the boundaries by using Imogen's multigrid summation algorithm. The linear
% solver is the biconjugate gradient method commonly used in gravitational hydrodynamics.
%
%>< run				Data manager                                            ImogenManager
%>< mass			Mass density                                            FluidArray
%<< phi 			Gravitational potential                                 GravityArray
        
    %--- Calculate BCs & Add Mass ---%
    %       Calculate the boundary conditions and add the mass for the right-hand side of the 
    %       Poisson's equation.
    if run.time.iteration < 4; tic; end

    bcsAndMass = calculateGravityEdge_GPU(mass, run.DGRID, run.gravity.mirrorZ, run.gravity.bconditionSource);
    if run.gravity.constant ~= 1
        bcsAndMass = run.gravity.constant*bcsAndMass;
    end

    if run.time.iteration < 4; t4bc = toc; end

    %--- Execute Linear Solver ---%
    %        phi0 is normally the potential from last step; This saves quite a bit of time.
    bcsAndMass = GPUdouble(bcsAndMass);
    bcsAndMass = L2_PolynomialPrecondition(bcsAndMass, size(mass.array), 0);

    phi0 = GPUdouble(phi0);

    [phi, flag, relres, iter] = bicgstab_GPU(@(x) findLaplacianTimesRHS(x, mass.gridSize, run.DGRID{1}), ...
                    bcsAndMass, run.gravity.tolerance, run.gravity.iterMax, ...
                    reshape(phi0, [prod(mass.gridSize) 1]) );

    if run.time.iteration < 4;
        fprintf('Phi solver (imogen step %i): BCs %.3gs, bicgstab_GPU in %.3gs/%.3g iter w/relres %6.6g\n', run.time.iteration, t4bc, toc-t4bc, iter, relres);
    end

    %--- Warn of Problems with Solver ---%
    if (flag)
        run.gravity.info = [run.gravity.info sprintf(['\nERROR - Gravity BiCgStab: ' ...
										'[Code Iteration %g] [Flag %g] ' ...
                                      '[Residual %g] [BiCgStab Iteration: %g of %g]'],...
                                      run.time.iteration, flag, relres, iter, run.gravity.iterMax)];
    end
    
    %--- Convert potential results back to domain-shaped array ---%
    phi = double(reshape(phi, mass.gridSize)); 

end

%--- This implements the 6th order HOC stencil of Spotz & Carey '95 w/o a super-huge matrix ---%
% This consumes most of the time spent by bicgstab_GPU
function Mx = findLaplacianTimesRHS(x, dims, h)

%--- These blocks of memory are used to store shifted-off edges in shiftAccumDaxpy ---%
persistent shiftStackX;
persistent shiftStackY;
persistent shiftStackZ;
persistent vSize;

if isempty(vSize) || (vSize ~= numel(x))
    N = max([dims(1)*dims(2) dims(2)*dims(3) dims(1)*dims(3) ]);
    shiftStackX = GPUdouble(zeros([N 1]));
    shiftStackY = GPUdouble(zeros([dims(1) 1 dims(3)]));
    shiftStackZ = GPUdouble(zeros([N 1]));
end

% HOC6 Laplacian prefactors
prefact = [-24 2 1 0] ./ (6*h*h);

x = reshape(x, dims);

Mx = GPUdouble(); setReal(Mx); setSize(Mx, dims); GPUallocVector(Mx);
GPUzeros(Mx);

%--- Perform circular-shift Laplacian operation ---%
for YSH = [-1 0 1];

    switch YSH
        case -1; assign(0, shiftStackY, x, [1 1 dims(1)], 1,   [1 1 dims(3)]); assign(1, x, 0, [1 1 dims(1)], 1  , [1 1 dims(3)]);
        case 1; assign(0, shiftStackY, x, [1 1 dims(1)], END, [1 1 dims(3)]); assign(1, x, 0, [1 1 dims(1)], END, [1 1 dims(3)]);
    end

    for XSH = [-1 0 1]; for ZSH = [-1 0 1];
        cenDistSqr = XSH^2 + YSH^2 + ZSH^2;
        if cenDistSqr == 3; continue; end

        shiftAccumDaxpy(x          , Mx, shiftStackX, shiftStackZ, [XSH YSH ZSH], prefact(cenDistSqr+1));
    end; end;

    switch YSH
        case -1; assign(1, x, shiftStackY, [1 1 dims(1)], 1  , [1 1 dims(3)]); assign(1, shiftStackY, 0, [1 1 dims(1)], 1  , [1 1 dims(3)]);
        case  1; assign(1, x, shiftStackY, [1 1 dims(1)], END, [1 1 dims(3)]); assign(1, shiftStackY, 0, [1 1 dims(1)], 1  , [1 1 dims(3)]);
    end

end

% Apply preconditioner (also reshapes)
Mx = L2_PolynomialPrecondition(Mx, dims, 1);
%Mx = reshape(Mx, [prod(dims) 1]);

end

function Ax = L2_PolynomialPrecondition(x, dims, noIniReshape)
% Polynomial preconditioner for Laplacians

persistent shiftStackX;
persistent shiftStackY;
persistent shiftStackZ;
persistent vSize;

if isempty(vSize) || (vSize ~= numel(x))
    N = max([dims(1)*dims(2) dims(2)*dims(3) dims(1)*dims(3) ]);
    shiftStackX = GPUdouble(zeros([N 1]));
    shiftStackY = GPUdouble(zeros([dims(1) 1 dims(3)]));
    shiftStackZ = GPUdouble(zeros([N 1]));
end

% Polynomials with 1-6 terms are effective. Polynomials with more do not help sufficiently.
numPolyTerms = 6;
polycoeffs = ...
[ 1      0      0      0      0      0      0      0      0      0      0      0; ...
  1.1660 0.8333 0      0      0      0      0      0      0      0      0      0; ...
  1.0930 1.5630 1.0930 0      0      0      0      0      0      0      0      0; ...
  0.9250 1.2250 2.2750 1.575  0      0      0      0      0      0      0      0; ...
  0.6789 0.0694 0.4155 2.4365 2.0969 0      0      0      0      0      0      0; ...
  1      0.5046 0.7589 0.9939 2.9454 2.9431 0      0      0      0      0      0; ...
  1      0.3976 0.6940 1.2421 1.5816 3.0061 2.2642 0      0      0      0      0; ...
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

Bn_x = GPUdouble(); setReal(Bn_x); setSize(Bn_x, dims); GPUallocVector(Bn_x);
cublasScopy(2*prod(dims), getPtr(Ax), 1, getPtr(Bn_x), 1);

accum = GPUdouble(); setReal(accum); setSize(accum, dims); GPUallocVector(accum);

%--- Iterate up to requested polynomial order ---%
for ord = 2:numPolyTerms;
    GPUzeros(accum);

    %--- Perform shift/accumulate operation ---%
    for XSH = [-1 0 1]; for YSH = [-1 0 1]; for ZSH = [-1 0 1];
        if (XSH^2 + YSH^2 + ZSH^2) == 1;
        % For basic case only act when shifting by exactly 1 (face-contact only)
            switch YSH
                case -1;
                    assign(0, shiftStackY, Bn_x, [1 1 dims(1)], 1,   [1 1 dims(3)]);
                    assign(1, Bn_x, 0, [1 1 dims(1)], 1  , [1 1 dims(3)]);
                    shiftAccumDaxpy(Bn_x, accum, shiftStackX, shiftStackZ, [XSH YSH ZSH], 1/6);
                    assign(1, Bn_x, shiftStackY, [1 1 dims(1)], 1  , [1 1 dims(3)]);
                case  0; shiftAccumDaxpy(Bn_x  , accum, shiftStackX, shiftStackZ, [XSH YSH ZSH], 1/6);
                case  1;
                    assign(0, shiftStackY, Bn_x, [1 1 dims(1)], END, [1 1 dims(3)]);
                    assign(1, Bn_x, 0, [1 1 dims(1)], END, [1 1 dims(3)]);
                    shiftAccumDaxpy(Bn_x, accum, shiftStackX, shiftStackZ, [XSH YSH ZSH], 1/6);
                    assign(1, Bn_x, shiftStackY, [1 1 dims(1)], END, [1 1 dims(3)]);
            end
        end
    end; end; end
    %--- Store next prefactors & accumulate next term in series ---%
    cublasScopy(2*Npoints, getPtr(accum), 1, getPtr(Bn_x), 1);
    shiftAccumDaxpy(accum, Ax, shiftStackX, shiftStackZ, [0 0 0], polycoeffs(numPolyTerms,ord));
end

Ax = reshape(Ax, [prod(dims) 1]);
end
