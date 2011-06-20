function phi = bicgstabPotentialSolver_GPU(run, mass, gridsize, iamrecursed)
% This routine solves for gravitational potential using a linear solver on the interior and 
% finding potential at the boundaries by using a multigrid summation algorithm. The linear
% solver is the biconjugate gradient method commonly used in gravitational hydrodynamics.
%
%>< run				Data manager                                            ImogenManager
%>< mass			Mass density                                            FluidArray
%<< phi 			Gravitational potential                                 GravityArray
        
    %--- Calculate BCs & Add Mass ---%
    %       Calculate the boundary conditions and add the mass for the right-hand side of the 
    %       Poisson's equation.
    if run.time.iteration < 4; tic; end

    % Recursively coarsen in order to propagate low-mode data faster
    if prod(gridsize) > 64^3 % If we > 64^3 on a side
        mprime = interpolateGPUvar(GPUdouble(mass), -2);
        philo = GPUdouble(bicgstabPotentialSolver_GPU(run, mprime, gridsize/2,1))/2;

        philo = reshape(philo, gridsize/2);
        phi0 = interpolateGPUvar(philo, 2);
        phi0 = reshape(phi0, [numel(phi0) 1]);
	else
        phi0 = zeros([numel(mass) 1], GPUdouble);
    end

    bcsAndMass = calculateGravityEdge_GPU(mass, run.DGRID, run.gravity.mirrorZ);
    if run.gravity.constant ~= 1
        bcsAndMass = run.gravity.constant*bcsAndMass;
    end

    if run.time.iteration < 4; t4bc = toc; end

    %--- Execute Linear Solver ---%
    %        phi0 is normally the potential from last step; This saves quite a bit of time.
    bcsAndMass = polyPreconditionL2_GPU(bcsAndMass, size(mass), 0);

    [phi, flag, relres, iter] = bicgstab_GPU(@(x) findLaplacianTimesRHS(x, gridsize, run.DGRID{1}), ...
                    bcsAndMass, run.gravity.tolerance, run.gravity.iterMax, phi0);

    if (run.time.iteration < 4)
	if nargin == 4; fprintf('    '); end
        fprintf('Phi solver (imogen step %i): BCs %.3gs, bicgstab_GPU in %.3gs/%.3g iter w/relres %6.6g\n', run.time.iteration, t4bc, toc-t4bc, iter, relres);
    end

    %--- Warn of Problems with Solver ---%
    if (flag)
        run.gravity.info = [run.gravity.info sprintf(['\nERROR - Gravity BiCgStab: ' ...
										'[Code Iteration %g] [Flag %g] ' ...
                                      '[Residual %g] [BiCgStab Iteration: %g of %g]'],...
                                      run.time.iteration, flag, relres, iter, run.gravity.iterMax)];
    end
    
    %--- Convert potential results back to domain-shaped array and cast back to double ---%
    phi = reshape(phi, gridsize); 
    if iamrecursed == 0; phi = double(phi); end

end

%--- This implements the 4th order HOC stencil of Spotz & Carey '95 in a matrix-free manner ---%
% This consumes most of the time spent by bicgstab_GPU
function Mx = findLaplacianTimesRHS(x, dims, h)

% HOC4 Laplacian prefactors
prefact = [-24 2 1 0] ./ (6*h*h);

% Reshape x and allocate array to store returned operation
x = reshape(x, dims);
Mx = GPUdouble(); setReal(Mx); setSize(Mx, dims); GPUallocVector(Mx);

% Call operator with HOC Laplacian weights and precondition it
symmetricLinearOperator(x, Mx, prefact);
Mx = polyPreconditionL2_GPU(Mx, dims, 1);
%Mx = reshape(Mx, [numel(Mx) 1]);

end

