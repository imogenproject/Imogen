classdef GravityManager < handle
% This is the management class for the potential solver. This is a singleton class to be accessed 
% using the getInstance() method and not instantiated directly. Currently the gravitational code is
% setup for a gravitational constant, G, of one.

    
%===================================================================================================
    properties (Constant = true, Transient = true) %                         C O N S T A N T     [P]
    end%CONSTANT
    
%===================================================================================================
    properties (SetAccess = public, GetAccess = public, Transient = true) %         P U B L I C  [P]
    ACTIVE;             % Specifies that gravity spolver state                          logical 
    info;               % Stores gravity solver report information                      str        

    laplacianMatrix;   % 6th order discrete Laplacian for gravity solver                sparse
    lowerConditioner;  % Incomplete LU factorizations to precondition the gravity       sparse
    upperConditioner;  % Solver for rapid solution                                      sparse

    constant;           % Gravitational scaling constant. (defaulted to 1)              double
    iterMax;            % Max number of iterations before gravity solver stops          double
    tolerance;          % Escape tolerance for the iterative solver                     double

    fixedPotential;     % Adds a predefined fixed potential to any solved-for potential sparse

    solve;              % Function handle to the potential solver                       handle
    initialize;         % Initialize function handle for gravitation solver             handle

    bconditionSource;   % Determines use of full or interpolated boundary conditions    string

    mirrorZ            % If true creates BCs with mass mirrored across lower XY plane  bool [false]

    TYPE;               % Gravity solver type enumeration                               str
    end%PUBLIC
    
%===================================================================================================
    properties (SetAccess = public, GetAccess = private) %                         P R I V A T E [P]
        pMatrixActive;     % Specifies sparse matrix active/inactive                    logical
        parent;            % Manager parent                                             ImogenManger
    end %PROPERTIES 
    
%===================================================================================================
    methods %                                                                     G E T / S E T  [M]
    end%GET/SET    
    
%===================================================================================================
    methods (Access = public) %                                                     P U B L I C  [M]        
        
%___________________________________________________________________________________________________ createSparseMatrix
% Builds a Laplacian and a preconditioner for linear solver methods. First creates the 6th order
% approximate Laplacian matrix and then performs a block incomplete LU factorization to create a
% preconditioner. The block size is chosen as Nx*Ny/2 as this has been found to be a good compromise
% between preconditioner size and facilitating rapid convergence.
function createSparseMatrix(obj, grid, dgrid)
            if (~obj.ACTIVE || ~obj.pMatrixActive); return; end

            fprintf('Generating subset of Laplacian matrix to build preconditioner...\n');
            
            blockSize = grid(1)*grid(2)/2;
            fprintf('Generating block ILU preconditioner, block size %i...\n', blockSize);
            [obj.lowerConditioner obj.upperConditioner] = ...
                poissonBlockILU(obj.laplacianMatrix, .05, blockSize, [prod(grid) prod(grid)]);

        end
%___________________________________________________________________________________________________ setSolver
% Attaches the correct solver function to the solve handle property as specified by the input type.
        function setSolver(obj, type)
            obj.TYPE = type;
            
            switch type
                %-----------------------------------------------------------------------------------
                case ENUM.GRAV_SOLVER_EMPTY
                    obj.solve           = @emptyPotentialSolver;
                    obj.initialize      = @emptyPotentialSolverIni;
                    obj.pMatrixActive   = false;
                %-----------------------------------------------------------------------------------
                case ENUM.GRAV_SOLVER_BICONJ
                    obj.solve           = @bicgstabPotentialSolver;
                    obj.initialize      = @bicgstabPotentialSolverIni;
                    obj.pMatrixActive   = true;
                case ENUM.GRAV_SOLVER_GPU
                    obj.solve           = @bicgstabPotentialSolver_GPU;
                    obj.initialize      = @bicgstabPotentialSolverIni_GPU;
                    obj.pMatrixActive   = false;
               %-----------------------------------------------------------------------------------
                case ENUM.GRAV_SOLVER_MULTIGRID
                    obj.solve           = @multigridPotentialSolver;
                    obj.initialize      = @multigridPotentialSolverIni;
                    obj.pMatrixActive   = false;
                %-----------------------------------------------------------------------------------
                otherwise
                    obj.type            = ENUM.GRAV_SOLVER_EMPTY;
                    obj.solve           = @emptyPotentialSolver;
                    obj.initialize      = @emptyPotentialSolverIni;
                    obj.pMatrixActive   = false;
            end
                
            
        end
        
%___________________________________________________________________________________________________ solvePotential
% Actual method call for finding the gravitational potential for a given mass distribution. Solver
% has an initial abort statement that exits if the gravitational solver is not active for a run.
        function solvePotential(obj, run, mass, grav)
            if ~obj.ACTIVE; return; end
            
            if (run.time.iteration == 0) && strcmp(obj.TYPE, ENUM.GRAV_SOLVER_EMPTY)
                grav.array = obj.fixedPotential;
            end

            if ~strcmp(obj.TYPE, ENUM.GRAV_SOLVER_EMPTY)
                grav.array = obj.fixedPotential + obj.solve(run, mass.array, mass.gridSize, 0);
            end


        end

%        function mainInitialize(obj, mass, grav)
            % 
 %       end

    end%PUBLIC
     
%===================================================================================================    
    methods (Access = private) %                                                P R I V A T E    [M]
        
%___________________________________________________________________________________________________ GravityManager
% Creates a new GravityManager instance and intializes it with default settings.
        function obj = GravityManager() 
            obj.setSolver( ENUM.GRAV_SOLVER_EMPTY );
            obj.ACTIVE      = false;
            obj.initialize  = @grav_ini_nonGravitational;

            obj.tolerance   = 1e-10;
            obj.iterMax     = 100;
            obj.constant    = 1;
            obj.fixedPotential = 0;

            obj.bconditionSource = ENUM.GRAV_BCSOURCE_FULL;
        end
        
    end%PROTECTED
        
%===================================================================================================    
    methods (Static = true) %                                                      S T A T I C    [M]
        
%___________________________________________________________________________________________________ getInstance
% Accesses the singleton instance of the GravityManager class, or creates one if none have
% been initialized yet.
        function singleObj = getInstance()
            persistent instance;
            if isempty(instance) || ~isvalid(instance) 
                instance = GravityManager();
            end
            singleObj = instance;
        end
        
    end%STATIC
    
end
