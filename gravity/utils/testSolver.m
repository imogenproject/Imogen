% This test script produces resolution/error outputs to test the potential solver's accuracy.
run.time.iteration = 0;
run.DGRID={1,1,1};
run.gravity.bconditionSource='full';
run.gravity.tolerance = 1e-10;
run.gravity.iterMax = 250;
run.gravity.constant = 1;
run.gravity.info='eh';
run.gravity.mirrorZ = 0;

%R = [16 32 64 128 256];
R = [128];
enorm = zeros(size(R));

for a = 1:numel(R);
    u = R(a);

    [X Y Z] = ndgrid(1:u,1:u,1:u);
    X = X-u/2;
    Y = Y-u/2;
    Z = Z-u/2;

    rad = 2*sqrt(X.^2+Y.^2+Z.^2)/u;

    mass.array = 1 - 2*rad + rad.^2;
    mass.array(rad > 1) = eps;

    mass.gridSize = size(mass.array);

    run.DGRID={2/u,2/u,2/u};

    phi = double(bicgstabPotentialSolver_GPU(run, mass));

    ALPHA=2*pi/3; BETA=-ALPHA; GAMMA=pi/5;

    phiAnalytic = ALPHA*rad.^2 + BETA*rad.^3 + GAMMA*rad.^4 - pi/3; % Constant to be in Coulomb gauge
    phiAnalytic(rad > 1) = -2*pi./(15*rad(rad > 1));
 
    errormat = phi - phiAnalytic;
    errormat(isnan(errormat)) = 0; % Avoid any singularities
    errormat(isinf(errormat)) = 0;

    enorm(a) = sqrt(sum(errormat(:).^2)/numel(mass.array));
end

