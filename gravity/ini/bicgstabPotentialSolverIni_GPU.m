function bicgstabPotentialSolverIni_GPU(run, mass, grav)
% This function initializes the bicgstab solver; It prepares the coefficient matrix and
% the preconditioner
%
%>< run         Data manager                                                    ImogenManager
%>< mass        Mass density                                                    FluidArray
%>< grav        Gravitational potential                                         GravityArray

%run.gravity.createSparseMatrix(mass.gridSize, [run.DGRID{1} run.DGRID{2} run.DGRID{3}]);

% FIXME: This needs to start up GPUmat.
% FIXME: GPUmat needs to be start-able without user interaction

grav.array = zeros(mass.gridSize);

end
