classdef Radiation < handle
% Contains the functionality for handling radiation based sources and sinks.
    
%===================================================================================================
    properties (Constant = true, Transient = true) %                            C O N S T A N T  [P]
    end%CONSTANT
    
%===================================================================================================
    properties (SetAccess = public, GetAccess = public, Transient = true) %         P U B L I C  [P]
        type;           % Enumerated type of radiation model to use.                string
        strength;       % Strength coefficient for the radiation calculation.       double
        exponent;       % Radiation exponent.                                       double
        initialMaximum; % Initial maximum radiation value used to calculate the     double
                        %   strength coefficient paramter.                          
        solve;          % Handle to raditaion function for simulation.              @func
    end%PUBLIC
    
%===================================================================================================
    properties (SetAccess = public, GetAccess = private) %                        P R I V A T E  [P]
    end %PRIVATE
    
    
    
    
    
    
%===================================================================================================
    methods %                                                                     G E T / S E T  [M]
        
%___________________________________________________________________________________________________ Radiation
% Creates a new Radiation instance.
        function obj = Radiation() 
            obj.strength = 1;
            obj.exponent = 0.5;
            obj.type     = ENUM.RADIATION_NONE;
        end
        
    end%GET/SET
    
%===================================================================================================
    methods (Access = public) %                                                     P U B L I C  [M]

%___________________________________________________________________________________________________ preliminary
        function preliminary(obj)
            switch (obj.type)
                %-----------------------------------------------------------------------------------
                case ENUM.RADIATION_NONE
                    obj.solve = @obj.noRadiationSolver;
                %-----------------------------------------------------------------------------------
                case ENUM.RADIATION_OPTICALLY_THIN
                    obj.solve = @obj.opticallyThinSolver;
                %-----------------------------------------------------------------------------------
                otherwise
                    obj.type  = ENUM.RADIATION_NONE;
                    obj.solve = @obj.noRadiationSolver;
            end
        end

%___________________________________________________________________________________________________ 
% Initialize the radiation solver parameters.
        function initialize(obj, run, mass, mom, ener, mag)
            if strcmp(obj.type, ENUM.RADIATION_NONE)
                obj.strength = 0;
                return
            end
                                
            unscaledRadiation = obj.solve(run, mass, mom, ener, mag);
            
            kineticEnergy     = 0.5*(mom(1).array .* mom(1).array + mom(2).array .* mom(2).array ...
                                    + mom(3).array .* mom(3).array) ./ mass.array;
            magneticEnergy    = 0.5*(mag(1).array .* mag(1).array + mag(2).array .* mag(2).array ...
                                    + mag(3).array .* mag(3).array);
            
            obj.strength      = (obj.initialMaximum/100)* ...
                   minFinderND( (ener.array - kineticEnergy - magneticEnergy) ./ unscaledRadiation );
        end
        
%___________________________________________________________________________________________________ noRadiationSolver
% Empty solver for non-radiating cases. 
        function result = noRadiationSolver(obj, run, mass, mom, ener, mag)
            result = 0;
        end
        
%___________________________________________________________________________________________________ opticallyThinSolver
% Solver for free radiation.
        function result = opticallyThinSolver(obj, run, mass, mom, ener, mag)
            gasPressure = pressure(ENUM.PRESSURE_GAS, run, mass, mom, ener, mag);
            result      = obj.strength*mass.array.^(2 - obj.exponent) .* gasPressure.^obj.exponent;
        end
        
    end%PUBLIC
    
    
%===================================================================================================    
    methods (Access = private) %                                                P R I V A T E    [M]        
    end%PROTECTED
        
%===================================================================================================    
    methods (Static = true) %                                                     S T A T I C    [M]        
    end%STATIC
    
end%CLASS
        
