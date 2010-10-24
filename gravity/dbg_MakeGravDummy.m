% Generates an environment where an arbitrary mass distribution can be tested using the Imogen solvers

testgrid = input('Make sure density is in mass.array. Abort and do so if this is not the case: ');

mass.gridSize = size(mass.array);
if numel(mass.gridSize) == 2
    mass.gridSize(3) = 1;
end

run.gravity.bconditionSource = 'full';
run.gravity.tolerance = 1e-10;
run.gravity.iterMax = 100;
run.gravity.constant = 1;
run.gravity.info = 'blah';

run.time.iteration = 0;

run.DGRID = {1,1,1};

%run.gravity.LAPLACIAN_MATRIX = createLaplacianMatrix(mass.gridSize);
%bsize = mass.gridSize(1)*mass.gridSize(2)/2;

%[run.gravity.LOWER_CONDITIONER run.gravity.UPPER_CONDITIONER] = poissonBlockILU(run.gravity.LAPLACIAN_MATRIX, .005, bsize, [prod(mass.gridSize) prod(mass.gridSize)]);

clear bsize;
clear testgrid;

