classdef FluidManager < handle
% Manager class for fluid fluxing operations and settings. This is a singleton class to be accessed 
% using the getInstance() method and not instantiated directly.
    
%===================================================================================================
    properties (Constant = true, Transient = true) %                            C O N S T A N T  [P]
    end%CONSTANT
    
%===================================================================================================
    properties (SetAccess = public, GetAccess = public, Transient = true) %         P U B L I C  [P]
        MINMASS;                  % Minimum allowed mass value.                         double
        MASS_THRESHOLD;           % Threshold above which gravity solver functions.     double
        ACTIVE;                   % Specifies the fluid fluxing solver state.           logical 
        freezeSpd;                % Freeze speed array.                                 FluxArray
        freezeSpdTVD;             % Freeze speed array for 2nd order fluxing.           FluxArray
        thresholds;               % Threshold values for gravitational fluxing.         struct
        viscosity;                % Artificial viscosity object.                        ArtificialViscosity
        radiation;                % Radiation object.                                   Radiation
        limiter;                  % Flux limiters to use for each flux direction.       cell(3)
    end%PUBLIC
    
%===================================================================================================
    properties (SetAccess = public, GetAccess = private) %                        P R I V A T E  [P]
        parent;            % Parent manager                                         ImogenManager
    end %PRIVATE
    
    
    
    
    
    
%===================================================================================================
    methods %                                                                     G E T / S E T  [M]
    end%GET/SET
    
    
%===================================================================================================
    methods (Access = public) %                                                     P U B L I C  [M]

%___________________________________________________________________________________________________ preliminary
        function preliminary(obj)
            obj.viscosity.preliminary();
            obj.radiation.preliminary();
        end
        
%___________________________________________________________________________________________________ createFreezeArray
% Initializes the freeze speed arrays for relaxed fluxing of the fluid variables.
        function createFreezeArray(obj)
            obj.freezeSpd     = FluxArray.empty(3,0);
            obj.freezeSpdTVD = FluxArray.empty(3,0);
            for i=1:3
                obj.freezeSpd(i)    = FluxArray(ENUM.VECTOR(i), FluxArray.FREEZE, ...
                                                obj.parent);
                obj.freezeSpdTVD(i) = FluxArray(ENUM.VECTOR(i), FluxArray.FREEZE, ...
                                                obj.parent);
            end
        end
        
%___________________________________________________________________________________________________ setFluxLimiters
% Initializes the flux limiters array based on the enumerated input.
        function setFluxLimiters(obj, limiters)
            obj.limiter = cell(1,3);
            fields      = {'x', 'y', 'z'};
            for i=1:3
                if isfield(limiters, fields{i});
                    obj.limiter{i} = obj.parseFluxLimiterEnum(limiters.(fields{i}));
                else
                    obj.limiter{i} = @vanleerLimiter;
                end
            end
        end
        
    end%PUBLIC
    
%===================================================================================================    
    methods (Access = private) %                                                P R I V A T E    [M]
        
%___________________________________________________________________________________________________ parseFluxLimiterEnum
        function result = parseFluxLimiterEnum(obj, limiterEnum)
            switch limiterEnum
                case FluxLimiterEnum.VAN_LEER,      result = @vanleerLimiter;
                case FluxLimiterEnum.SUPERBEE,      result = @superbeeLimiter;
                case FluxLimiterEnum.MINMOD,        result = @minmodLimiter;
                otherwise
                    error(['Imogen:UnknownType: Unknown flux limiter type: ', limiterEnum]);
            end
        end
    
%___________________________________________________________________________________________________ FluidManager
% Creates a new FluidManager instance.
        function obj = FluidManager() 
            obj.ACTIVE                   = true;
            obj.MASS_THRESHOLD           = 0;
            obj.viscosity                = ArtificialViscosity();
            obj.radiation                = Radiation();
        end
        
    end%PROTECTED
        
%===================================================================================================    
    methods (Static = true) %                                                     S T A T I C    [M]
        
%___________________________________________________________________________________________________ getInstance
% Accesses the singleton instance of the FluidManager class, or creates one if none have
% been initialized yet.
        function singleObj = getInstance()
            persistent instance;
            if isempty(instance) || ~isvalid(instance) 
                instance = FluidManager(); 
            end
            singleObj = instance;
        end
        
    end%STATIC
    
end%CLASS
        
