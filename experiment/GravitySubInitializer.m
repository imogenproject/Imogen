classdef GravitySubInitializer < handle
% Handles all of the various initialization properties for the gravitational functionality within
% Imogen.
        
%===================================================================================================
        properties (Constant = true, Transient = true) %                 C O N S T A N T         [P]
    end%CONSTANT
        
%===================================================================================================
    properties (SetAccess = public, GetAccess = public) %                           P U B L I C  [P]
        fixedPotential;     % Contains a predefined, unchanging potential          double
        constant;           % Gravitational constant (G).                          double
        solver;             % Enumeration specifying which solver to use.          string
        iterMax;            % Maximum number of iterations by a linear technique   int
        tolerance;          % Tolerance for linear solver techniques               double
        bconditionSource;   % String: Determines function used to find BCs         string
        mirrorZ;            % If 1, mirrors gravity across the lower Z plane       true/false

        laplacianMatrix;    % These aren't really used anymore
        lowerConditioner;
        upperConditioner;

    end %PUBLIC

%===================================================================================================
    properties (SetAccess = protected, GetAccess = protected) %                P R O T E C T E D [P]
        pIniPointPotentials; % A cell array of point potentials to add to the simulation
        pNumIniPointPots;

    end %PROTECTED


        
        
%===================================================================================================
    methods %                                                                     G E T / S E T  [M]
        
%___________________________________________________________________________________________________ GravitySubInitializer
% Creates a new GravitySubInitializer object and sets the default settings.
        function obj = GravitySubInitializer()
            obj.constant         = 1;
            obj.solver           = ENUM.GRAV_SOLVER_EMPTY;
            obj.fixedPotential   = 0;
            obj.iterMax          = 150;
            obj.tolerance        = 1e-10;
            obj.bconditionSource = ENUM.GRAV_BCSOURCE_FULL;
            obj.mirrorZ          = 0;
        end
        
        end%GET/SET
        
%===================================================================================================
    methods (Access = public) %                                                     P U B L I C  [M]
        
    end%PUBLIC
        
%===================================================================================================        
        methods (Access = protected) %                                      P R O T E C T E D    [M]
        end%PROTECTED
                
%===================================================================================================        
        methods (Static = true) %                                                                                                          S T A T I C    [M]
        end%PROTECTED
        
end%CLASS
