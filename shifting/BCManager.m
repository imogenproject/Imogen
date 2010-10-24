classdef BCManager < handle
% This is the class responsible for managing artificial boundary conditions for Imogen runs. This is
% a  singleton class to be accessed using the getInstance() method and not instantiated directly.


%===================================================================================================
    properties (Constant = true, Transient = true) %                            C O N S T A N T  [P]
    end%CONSTANTS
    
%===================================================================================================
    properties (SetAccess = public, GetAccess = public, Transient = true) %         P U B L I C  [P]
        infinity;   % Number of cells to define infinity for bc interpolation.        int
        modes;      % Boundary condition modes                                        struct
    end%PUBLIC
    
%===================================================================================================
    properties (SetAccess = public, GetAccess = private) %                        P R I V A T E  [P]
        parent;     % Parent manager                                            ImogenManager
    end %PRIVATE

    
    
    
    
    
    
%===================================================================================================
    methods %                                                                     G E T / S E T  [M]
    end%GET/SET
    
    
%===================================================================================================
    methods (Access = public) %                                                     P U B L I C  [M]
      
%___________________________________________________________________________________________________ attachBoundaryConditions
% Reads the boundary conditions data object (either a structure or array) from the initialization
% properties set in the BCManager object.
%>< arrayObj    ImogenArray object on which to attach boundary conditions.          ImogenArray
        function attachBoundaryConditions(obj, arrayObj)
            %--- Same bcs for all case ---%
            %       If no explicit bounardies conditions were set for the arrayObj, the general
            %       boundary conditions are applied.
            if ~isfield(obj.modes, arrayObj.id{1})
                obj.populateBCModes(arrayObj, obj.modes);
                return
            end   
               
            %--- Handle FluxArray cases ---%
            %       Flux boundary conditions can be explicitily set separately from the primary 
            %       variable arrays. Hence they need to be handled accordingly, which requires two
            %       levels of id verification 1) for the primary array to which they belong and
            %       2) for the kind of flux they represent.
            if isa(arrayObj,'FluxArray')
                
                %--- Flux specific (targeted) case ---%
                if isfield(obj.modes.(arrayObj.id{1}), arrayObj.id{2})    
                    obj.populateBCModes(arrayObj, obj.modes.(arrayObj.id{1}).(arrayObj.id{2}), ...
                                        obj.modes.(arrayObj.id{1}), obj.modes);

                %--- Global flux case ---%
                elseif isfield(obj.modes, arrayObj.id{2})
                    obj.populateBCModes(arrayObj, obj.modes.(arrayObj.id{2}), obj.modes);

                %--- Primary variable inherited case ---%
                else
                    obj.populateBCModes(arrayObj, obj.modes.(arrayObj.id{1}), obj.modes);
                end
                
                return
            end

            %--- Remaining cases ---%
            %       Handles the boundaries for the remaining array types both primary variables and
            %       other non-flux objects like freezing speed, grid velocity, etc.
            obj.populateBCModes(arrayObj, obj.modes.(arrayObj.id{1}), obj.modes);

        end
                
%___________________________________________________________________________________________________ toStruct
% Converts the pertinent BCManager data to a BC structure for saving or non-class use.
% # result        The BCManager data output as a structure.                                Struct
        function result = toStruct(obj)
            result.modes    = obj.modes;
            result.infinity = obj.infinity;
        end
        
    end%PUBLIC
    
%===================================================================================================
    methods (Access = private) %                                                P R I V A T E    [M]
    
%___________________________________________________________________________________________________ populateBCModes
% Responsible for actually initializeing the boundary conditions on a per boundary condition
% mode type. Once it has initialized the mode it applies the boundary condition to the 
% ImogenArray by setting the bcModes value and attaching the edgeshift routines.
%>> varargin    Variable length input from the readBoundaryConditions method.        cell
        function populateBCModes(obj, arrayObj, varargin)
            bcs                 = varargin;
            result.x            = [];  
            result.y            = []; 
            result.z            = [];
            fields              = {'x','y','z'};
            arrayObj.bcModes    = cell(2,3);
            
            for i=1:3
                for n=1:length(bcs)
                    if ( isfield(bcs{n},fields{i}) && isempty(result.(fields{i})) )
                        result.(fields{i}) = bcs{n}.(fields{i});
                        break
                    end
                end
            end
            
            %--- Check for empty boundary conditions ---%
            for i=1:3
                if isempty(result.(fields{i}))
                    warning('ImogenArray:BCModeINIError',[arrayObj.idAsString() ...
                        ' is missing boundary condition for ' fields{i} ...
                        ' direction. Using default CIRCULAR condition instead.']);
                    result.(fields{i}) = obj.modes{i};
                end
            end
            
            for i=1:3
                if iscell(result.(fields{i}))
                    secondIndex = length(result.(fields{i}));
                    arrayObj.bcModes{1,i} = result.(fields{i}){1};
                    arrayObj.bcModes{2,i} = result.(fields{i}){secondIndex};
                else
                    arrayObj.bcModes{1,i} = result.(fields{i});
                    arrayObj.bcModes{2,i} = result.(fields{i});
                end
            end

            arrayObj.edgeshifts  = cell(2,3);
            for i=1:3
                for n=1:2
                    switch arrayObj.bcModes{n,i}
                        case ENUM.BCMODE_CIRCULAR,    arrayObj.edgeshifts{n,i} = @circ_shift;
                        case ENUM.BCMODE_CONST,       arrayObj.edgeshifts{n,i} = @constant_shift;
                        case ENUM.BCMODE_FADE,        arrayObj.edgeshifts{n,i} = @fade_shift;
                        case ENUM.BCMODE_FLIP,        arrayObj.edgeshifts{n,i} = @flip_shift;
                        case ENUM.BCMODE_LINEAR,      arrayObj.edgeshifts{n,i} = @linear_shift;
                        case ENUM.BCMODE_MIRROR,      arrayObj.edgeshifts{n,i} = @mirror_shift;
                        case ENUM.BCMODE_TRANSPARENT, arrayObj.edgeshifts{n,i} = @transparent_shift;
                        case ENUM.BCMODE_WALL,        arrayObj.edgeshifts{n,i} = @wall_shift;
                        case ENUM.BCMODE_WORMHOLE,    arrayObj.edgeshifts{n,i} = @wormhole_shift;
                        case ENUM.BCMODE_ZERO,        arrayObj.edgeshifts{n,i} = @zero_shift;
                        otherwise
                            error(['Imogen:UnknownType: Unknown BC mode ', arrayObj.bcModes{n,i}]);
                    end
                end
            end
        end
        
%___________________________________________________________________________________________________ BCManager
% Creates a new BCManager instance.
        function obj = BCManager(); end
        
    end%PROTECTED
        
%===================================================================================================    
    methods (Static = true) %                                                      S T A T I C   [M]
        
%___________________________________________________________________________________________________ getInstance
% Accesses the singleton instance of the BCManager class, or creates one if none have
% been initialized yet.
        function singleObj = getInstance()
            persistent instance;
            if isempty(instance) || ~isvalid(instance) 
                instance = BCManager(); 
            end
            singleObj = instance;
        end
      
    end%STATIC
end%CLASS
