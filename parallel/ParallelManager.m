classdef ParallelManager < handle
% Manages the parallel settings and operations in the Imogen code. This is a singleton class to be 
% accessed using the getInstance() method and not instantiated directly.
    
%===================================================================================================
    properties (Constant = true, Transient = true) %                            C O N S T A N T  [P]
    end%CONSTANT
    
%===================================================================================================
    properties (SetAccess = public, GetAccess = public) %                            P U B L I C [P]
        ACTIVE;             % Specifies the ACTIVE state of parallelization.        logical
        poolSize;           % CPUs used for parallel distribution.                  int
        distribution;       % Specifies how to distribute parallel arrays.          codistributor
        parent;             % Parent manager.                                       ImogenManager
    end %PUBLIC
        
%===================================================================================================
    properties (SetAccess = private, GetAccess = private) %                        P R I V A T E [P]
        pRegisteredArrays;	% List of arrays to include in redistribution.          cell(?)
        pRegisteredIndex;   % Number of array registerd for redistribution.         int
        pDimensionPriority;  % Preferential order for codistribution dimensions.    int(3)
    end %PRIVATE
    
    
    
    
    
%===================================================================================================
    methods %                                                                     G E T / S E T  [M]        
    end%GET/SET
    
%===================================================================================================
    methods (Access = public) %                                                      P U B L I C [M]
        
%___________________________________________________________________________________________________ registerForRedistribution
% Adds an ImogenArray... object to the list of items to be redistributed
        function registerForRedistribution(obj,imoArray)
            obj.pRegisteredIndex                        = obj.pRegisteredIndex + 1;
            obj.pRegisteredArrays{obj.pRegisteredIndex} = imoArray;
        end
        
        
%___________________________________________________________________________________________________ redistributeArrays
% Changes the distribution scheme used by the arrays.
%    dim        dimension along which to splice the arrays
        function redistributeArrays(obj,dim)
            if (dim ~= obj.distribution.Dimension), return; end % Skip if no distribution conflict.
            
            obj.distribution = codistributor1d(obj.pDimensionPriority(dim));
            for i=1:obj.pRegisteredIndex
                obj.pRegisteredArrays{i}.redistribute(obj.distribution);
            end
        end
        
%___________________________________________________________________________________________________ preliminary
% Handles the preliminary initialization for parallel operation if running in parallel mode, or
% sets parameters to run properly in serial mode.
        function preliminary(obj)
            
            %--- Initialize parallel variables ---%
            if numlabs > 1 
                obj.ACTIVE = true;
                obj.parent.save.logPrint('Running in parallel mode (pool size = %g).\n',numlabs);
            else
                obj.parent.save.logPrint('Running in serial mode.\n');
                return;
            end
             
            %--- Dimensional Priorities ---%
            %       The most efficient codistribution of parallel arrays is found by maximizing the
            %       grid size along each prioritized direction.
            [maxGrid, maxIndex]     = max(obj.parent.gridSize);
            obj.pDimensionPriority  = maxIndex * ones(1,3);
            minGrid                 = min(obj.parent.gridSize);
            
            if sum(maxGrid == obj.parent.gridSize) > 1 % Case with two maximum dimensions.
                obj.pDimensionPriority(maxIndex) = ...
                                        find(obj.parent.gridSize == maxGrid && 1:3 ~= maxIndex, 1);
                                    
            elseif sum(minGrid == obj.parent.gridSize) > 1 % Case where two minimum dimensions.
                obj.pDimensionPriority(maxIndex) = ...
                                        find(obj.parent.gridSize == minGrid, 1, 'first');
                                    
            else % Case with dimensions of three distinct sizes.
                obj.pDimensionPriority(maxIndex) = ...
                                        find(obj.parent.gridSize > minGrid && 1:3 ~= maxIndex, 1);
            end

            %--- Set initial values ---%
            obj.pRegisteredIndex  = 0;
            obj.pRegisteredArrays = cell(1,20);
            obj.distribution      = codistributor1d(obj.pDimensionPriority(1));
        end
        
    end%PUBLIC
    
%===================================================================================================    
    methods (Access = private) %                                                P R I V A T E    [M]
        
%___________________________________________________________________________________________________ ParallelManager
        % Creates a new ParallelManager instance.
        function obj = ParallelManager() 
            obj.ACTIVE       = false;
            obj.distribution = codistributor1d();
        end
        
    end%PROTECTED
        
%===================================================================================================    
    methods (Static = true) %                                                      S T A T I C   [M]
        
%___________________________________________________________________________________________________ getInstance
% Accesses the singleton instance of the ParallelManager class, or creates one if none have been 
% initialized yet.
        function singleObj = getInstance()
            persistent instance;
            if isempty(instance) || ~isvalid(instance) 
                instance = ParallelManager(); 
            end
            singleObj = instance;
        end
        
    end%STATIC
    
end%CLASS
