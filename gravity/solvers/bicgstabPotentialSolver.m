function phi = bicgstabPotentialSolver(run, mass, phi0)
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

    bcsAndMass = calculateGravityEdge(mass, run.DGRID, run.gravity.mirrorZ, run.gravity.bconditionSource);
    if run.gravity.constant ~= 1
        bcsAndMass = run.gravity.constant*bcsAndMass;
    end

    if run.time.iteration < 4; t4bc = toc; end

    %--- Execute Linear Solver ---%
    %        phi0 is normally the potential from last step; This saves quite a bit of time.

    bcsAndMass = L2_PolynomialPrecondition(bcsAndMass, mass.gridSize, 0);

    [phi, flag, relres, iter] = bicgstab(@(x) findLaplacianTimesRHS(x, mass.gridSize, run.DGRID{1}), bcsAndMass, ...
                                    run.gravity.tolerance, run.gravity.iterMax, ...
                                    [], [], reshape(phi0, [numel(phi0) 1]) );

    if run.time.iteration < 4;
        fprintf('Phi solver stats (imogen iter %i): BCs in %.4gs, bicgstab in %.4gs/%i iter w/relres %g\n', run.time.iteration, t4bc, toc-t4bc, iter, relres);
    end

    %--- Warn of Problems with Solver ---%
    if (flag)
        run.gravity.info = [run.gravity.info sprintf(['\nERROR - Gravity BiCgStab: ' ...
										'[Code Iteration %g] [Flag %g] ' ...
                                      '[Residual %g] [BiCgStab Iteration: %g of %g]'],...
                                      run.time.iteration, flag, relres, iter, run.gravity.iterMax)];
    end
    
    %--- Convert potential results back to domain-shaped array ---%
    phi = reshape(phi, mass.gridSize); 

end

%--- This implements the 6th order HOC stencil of Spotz & Carey '95 w/o a super-huge matrix ---%
function Mx = findLaplacianTimesRHS(x, dims, h)
prefact = [-24 2 1 0] ./ (6*h*h);

x = reshape(x, dims);
Mx = zeros(size(x));

%--- Perform circular-shift Laplacian operation ---%
for XSH = [-1 0 1]; for YSH = [-1 0 1]; for ZSH = [-1 0 1];
   if (XSH^2 + YSH^2 + ZSH^2) == 3; continue; end
   q = circshift(x, [XSH YSH ZSH]);

    % Delete garbage in shifted array
    if XSH ==  1; q(1,:,:) = 0; end
    if XSH == -1; q(end,:,:) = 0; end

    if YSH ==  1; q(:,1,:) = 0; end;
    if YSH == -1; q(:,end,:) = 0; end;

    if ZSH ==  1; q(:,:,1) = 0; end;
    if ZSH == -1; q(:,:,end) = 0; end;

    switch (XSH^2 + YSH^2 + ZSH^2);
        case 0; Mx = Mx + prefact(1)*q; % center
        case 1; Mx = Mx + prefact(2)*q; % face-connected
        case 2; Mx = Mx + prefact(3)*q; % edge-connected
%        case 3; Mx = Mx + prefact(4)*q; % corner-connected
    end

end; end; end

%--- Reshape into a vector and return ---%
%Mx = reshape(Mx, [prod(dims) 1]);
Mx = L2_PolynomialPrecondition(Mx, dims, 1);

end


function Ax = L2_PolynomialPrecondition(x, dims, noIniReshape)
% Polynomial preconditioner for Laplacians
persistent shiftStackX;
persistent shiftStackY;
persistent shiftStackZ;

if isempty(shiftStackX)
    N = max([dims(1)*dims(2) dims(2)*dims(3) dims(1)*dims(3) ]);
    shiftStackX = zeros([N 1]);
    shiftStackY = zeros([dims(1) 1 dims(3)]);
    shiftStackZ = zeros([N 1]);
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
Ax  = polycoeffs(numPolyTerms,1) * x;
Bn_x = x;

%--- Iterate up to requested polynomial order ---%
for ord = 2:numPolyTerms;
    accum = zeros(size(x));

    %--- Perform shift/accumulate operation ---%
    for XSH = [-1 0 1]; for YSH = [-1 0 1]; for ZSH = [-1 0 1];
        if (XSH^2 + YSH^2 + ZSH^2) == 1;
            accum = accum + circshift(Bn_x, [XSH YSH ZSH]);
            
            if XSH ==  1; accum(1,:,:)   = accum(1,:,:) - Bn_x(end,:,:); end
            if XSH == -1; accum(end,:,:) = accum(end,:,:) - Bn_x(1,:,:); end
            if YSH ==  1; accum(:,1,:)   = accum(:,1,:) - Bn_x(:,end,:); end
            if YSH == -1; accum(:,end,:) = accum(:,end,:) - Bn_x(:,1,:); end
            if ZSH ==  1; accum(:,:,1)   = accum(:,:,1) - Bn_x(:,:,end); end
            if ZSH == -1; accum(:,:,end) = accum(:,:,end) - Bn_x(:,:,1); end
        
        end
    end; end; end

    %--- Store next prefactors & accumulate next term in series ---%
    Bn_x = accum / 6;
    Ax = Ax + polycoeffs(numPolyTerms,ord) * accum / 6;
end

Ax = reshape(Ax, [prod(dims) 1]);

end
