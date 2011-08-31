classdef ImogenArray < handle
% The base class for all Imogen array objects. Provides all of the common functionality shared by 
% the flux-source array objects.
    
%===================================================================================================
    properties (GetAccess = public, Constant = true, Transient = true) %        C O N S T A N T  [P]
    end%CONSTANT
    
%===================================================================================================
    properties (SetAccess = public, GetAccess = public) %                           P U B L I C  [P]
        component;      % Vector component of the array (0,1,2,3).                  int
        id;             % Identifier specifying the object data type.               cell(1,?)
        bcModes;        % Boundary condition types for each direction.              cell(2,3)
        bcInfinity;     % Number of cells to infinity for edges.                    int
        edgeStore;      % Stored edge values for shifting.                          Edges

        staticActive;   % Determines if statics should be applied.                  logical
        staticValues;   % Values for static array indices.                          double
        staticCoeffs;   % Coefficients used to determine how fast the array fades to the value double
        staticIndices;  % Indices of staticValues used by StaticArray.                int
        staticLinIndices;% The linear indexes into the array used to apply          int

        edgeshifts;     % Handles to shifting functions for each grid direction.    handle(2,3)
        isZero;         % Specifies that the array is statically zero.              logical

        indexGriddim;
    end %PROPERTIES
    
%===================================================================================================
    properties (SetAccess = protected, GetAccess = protected) %                P R O T E C T E D [P]
        pArray;         % Storage of the 3D array data.                             double(Nx,Ny,Nz)
        pShiftEnabled;  % Specifies if shifting is active for each dimension.       logical(3)
        pDistributed;   % Specifies if the object should be distributed.            bool
        pCompatibility; % Matlab version compability identifier.                    str
        pRunManager;    % Access to the ImogenManager singleton for the run.        ImogenManager
        pFades;         % Fade objects for array processing.                        cell(?)
        pFadesValue;    % Values to use for each of the fade objects.               double
        pUninitialized; % Specifies if obj has been initialized.                    bool
    end %PROPERTIES
   
%===================================================================================================
    properties (Dependent = true) %                                            D E P E N D E N T [P]
        gridSize;       % Size of the data array.                                   int(3)
        array;          % Storage of the 3D array data.                             double(Nx,Ny,Nz)
        idString;       % String form of the id cell array.                         str
        fades;          % The fade objects influencing the ImogenArray object.      cell{?}
    end %DEPENDENT

    
    
%===================================================================================================
    methods %                                                                     G E T / S E T  [M]

%___________________________________________________________________________________________________ ImogenArray
% Creates a new ImogenArray object according to the specified inputs.
%>> id          Identification information for the new object.                      cell/str
%>< run         Run manager object.                                                 ImogenManager
%>> statics     Static arrays and values structure.                                 struct
        function obj = ImogenArray(component, id, run, statics)
            if (nargin == 0 || isempty(id)); return; end
            if ~isa(id,'cell');    id = {id}; end
            
            obj.pCompatibility  = str2double(run.matlab.Version);
            obj.pShiftEnabled   = true(1,3);
            obj.component       = component;
            obj.id              = id;
            obj.bcInfinity      = run.bc.infinity;
            obj.pDistributed    = run.parallel.ACTIVE;
            obj.pRunManager     = run;
            obj.pFadesValue     = 0.995;
            obj.readStatics(statics);
            run.bc.attachBoundaryConditions(obj);
        end

%___________________________________________________________________________________________________ GS: fades
% Accesses the fades attached to the ImogenArray object.
        function result = get.fades(obj),   result = obj.pFades;        end
        
%___________________________________________________________________________________________________ GS: array
% Main property accessor to the data array for the ImogenArray object.
        function result = get.array(obj),   result = obj.pArray;        end
        
        function set.array(obj,value)
        % Sets the data array to the new value and cleans up faded and static cells.
            obj.pArray        = value;

            if ~isempty(obj.pFadesValue),       obj.applyFades();       end % Fade array.
            if obj.staticActive,                obj.applyStatics();     end % Enforce static values.
            if obj.pUninitialized,              obj.applyInitialize();  end % Initialize if needed.
        end
        
%___________________________________________________________________________________________________ GS: gridSize
% Specifies the grid size of the data array.
        function result = get.gridSize(obj)
            result = size(obj.pArray);
            if ( length(result) < 3), result = [result 1]; end
        end
        
        function set.gridSize(obj,value)
            warning('ImogenArray:GridSize',['The gridSize property is read-only and is ' ...
                    'currently %s and cannot be set to %s.'], mat2str(size(obj.pArray)), value);
        end

%___________________________________________________________________________________________________ GS: idString
        function result = get.idString(obj)
            result = obj.id{1};
            if length(obj.id) > 1
                for i=2:length(obj.id)
                    result = strcat(result,'_',obj.id{i});
                end
            end
        end
        
        function set.idString(obj,value)
            warning('ImogenArray:IdString',['The idString property is read-only and cannot be'...
                    'set to %s directly on %s.'],value, obj.id{:});
        end
                
        
    end%GET/SET
        
%===================================================================================================
    methods (Access = public) %                                                     P U B L I C  [M]
        
%___________________________________________________________________________________________________ distribute
% Distributes the serial data array into a data parallel environment according to the distributor
% argument.
%>< dist    Codistributor object that specifies the parallel distribution.          Codistributor
        function distribute(obj, dist)
            if (obj.pCompatibility > 7.8)
                if (~isreplicated(obj.pArray) || ~iscodistributed(obj.pArray))
                    obj.pArray = labBroadcast(1, obj.pArray);
                end
                obj.pArray  = codistributed(obj.pArray, dist);
            else % Handles the change in codistributed between r2009a and r2009b
                obj.pArray = codistributed(obj.pArray, dist, 'convert');
            end
        end
        
%___________________________________________________________________________________________________ redistribute
% Redistributes a parallel array according to the distributor argument.
%>< dist    Codistributor object the specifies the new distribution.                Codistributor
        function redistribute(obj, dist)
            obj.pArray = redistribute(obj.pArray,dist);
        end
        
%___________________________________________________________________________________________________ cleanup
% Cleans up the ImogenArray by emptying the data array, reducing memory requirements.
        function cleanup(obj),    obj.pArray =  [];    end
        
%___________________________________________________________________________________________________ shift
% Shifts the input array along the specified direction and by the specified number of cells
% according to the edgeshift function specified by the edgeshifts cell array.
% * DIRECT        The direction along which to shift (1,2,3)                  int
% * nCells        The number of cells to shift.                               int
% # result        The shifted array.                                          double    (Nx,Ny,Nz)
        function result = shift(obj, DIRECT, nCells)
            upperLowerIndex = 1 + (nCells > 0);

            if obj.pShiftEnabled(upperLowerIndex, DIRECT)
                result = obj.edgeshifts{upperLowerIndex, DIRECT}(obj.pArray, DIRECT, nCells, obj);
            else
                result = obj.pArray;
            end
            
            if (obj.staticActive); obj.applyStatics(); end
        end
        
%___________________________________________________________________________________________________ transparentEdge
% Returns the stored edge conditions used by transparent boundary conditions.
        function result = transparentEdge(obj,dim,upper)
            result = obj.edgeStore.getEdge(upper, dim, obj.pArray, obj.bcModes{1 + upper, dim});
        end

%___________________________________________________________________________________________________ idAsString
% Converts the ImogenArray object's identifier (id) property from a cell array to a string for 
% display or file writing purposes.
        function result = idAsString(obj)
            result = stringBuffer( ['[' class(obj) ']'] );
            for i=1:length(obj.id),     result.add(['.' obj.id{i}]); end
            result = result.string;
        end
        
%___________________________________________________________________________________________________ calculate5PtDerivative
% Calculates and returns the spatial derivative of the array.
% result = calculate5PtDerivative(obj,X,dGrid)
% X            Spatial direction in which to calculate the derivative.
% dGrid        Grid spacing in the X direction.
% result    Derivative of the potential array in the X direction.
        function result = calculate5PtDerivative(obj,X,dGrid)
            result = (    -obj.shift(X,2) + 8*obj.shift(X,1) ...
                        - 8*obj.shift(X,-1) + obj.shift(X,-2) ) ./ (12*dGrid);
        end    
        
%___________________________________________________________________________________________________ calculate2PtDerivative
% Calculates and returns the spatial derivative of the array.
% result = obj.calculate2PtDerivative(X,dGrid)
% X            Spatial direction in which to calculate the derivative.
% dGrid        Grid spacing in the X direction.
% result    Derivative of the potential array in the X direction.
        function result = calculate2PtDerivative(obj,X,dGrid)
            result = ( obj.shift(X,1) - obj.pArray ) ./ dGrid;
        end

%___________________________________________________________________________________________________ arrayIndexExchange
% Flips the array and all associated subarrays such that direction i and the x (stride-of-1) direction
% exchange places. Updates the array, all subarrays, and the static indices.
        function arrayIndexExchange(obj, toex, type)
            if numel(obj.indexGriddim) < 3; obj.indexGriddim = obj.gridSize; end

            if toex == 1; return; end

            l = [toex 2 3]; l(toex) = 1;
            obj.indexGriddim = obj.indexGriddim(l);
            
            if numel(obj.staticIndices) > 0
                ad = obj.indexGriddim;
                if toex == 2
                obj.staticIndices(:,2:4) = obj.staticIndices(:,[3 2 4]);
                obj.staticIndices(:,1)   = obj.staticIndices(:,2) + ...
                                          (obj.staticIndices(:,3)-1)*ad(1) + ...
                                          (obj.staticIndices(:,4)-1)*ad(1)*ad(2);
                end
                if toex == 3
                obj.staticIndices(:,2:4) = obj.staticIndices(:,[4 3 2]);
                obj.staticIndices(:,1)   = obj.staticIndices(:,2) + ...
                                          (obj.staticIndices(:,3)-1)*ad(1) + ...
                                          (obj.staticIndices(:,4)-1)*ad(1)*ad(2);
                end
            end

            if numel(obj.staticValues) > 0; obj.staticLinIndices = GPUdouble(obj.staticIndices(:,1)-1); end;
            if type == 1; obj.array = cudaArrayRotate(obj.array, toex); end

            %if strcmp(obj.id{1},'mag') && (numel(obj.id) == 1); obj.updateCellCentered(); end

        end

%___________________________________________________________________________________________________ applyStatics
% Applies the static conditions for the ImogenArray to the data array. This method is called during
% array assignment (set.array).
        function applyStatics(obj)
            if isa(obj.pArray,'double')
                if numel(obj.staticValues) > 0;
                    obj.pArray(obj.staticIndices(:,1)) = obj.pArray(obj.staticIndices(:,1)) + obj.staticCoeffs.*(obj.staticValues - obj.pArray(obj.staticIndices(:,1)));
                end
            else
                if numel(obj.staticValues) > 0; cudaStatics(obj.pArray, obj.staticLinIndices, obj.staticValues, obj.staticCoeffs, 8); end
            end
        end

        
    end%PUBLIC
    
%===================================================================================================
    methods (Access = protected) %                                          P R O T E C T E D    [M]

%___________________________________________________________________________________________________ applyInitialize
% Applies the initialization functions uninitialized ImogenArrays and sets the pUninitialized
% property to false to specify that the array has now been initialized. This method does little
% in the base class as it is meant to be extended for use by the InitializedArray subclass, which
% does not initialize during construction because the data is transient.
        function applyInitialize(obj),  obj.pUninitialized = false;     end
        
%___________________________________________________________________________________________________ applyFades
% Applies any fades in attached to the ImogenArray to the data array. This method is called during
% array assignment (set.array).
        function applyFades(obj)
            for i=1:length(obj.pFades)
                obj.pArray = obj.pFades{i}.fadeArray(obj.pArray, obj.pFadesValue);
            end
        end
        
%___________________________________________________________________________________________________ readFades
% Reads the fades stored in the ImogenManager object and applies the appropriate ones to this
% ImogenManager object. Note, fades cannot be applied to the FluxArray subclass because they result
% in unstable fluxing. So FluxArrays abort this property.
%>< run     Run manager object.                                                     ImogenManager
        function readFades(obj, run)
            if isa(obj,'FluxArray'); return; end
            index = 1;
            obj.pFadesValue = 0;
            %--- Find fade objects for this object ---%
            for i=1:length(run.fades)
                if any(strcmp(obj.id{1}, run.fades{i}.activeList))
                    obj.pFades{index} = run.fades{i};
                    index = index + 1;
                end
            end
            
        end
        
        
%___________________________________________________________________________________________________ initializeShiftingStates
% Determines which grid dimensions have extent, i.e. are greater than 1 cell, and activates shifting
% on only those dimensions to speed up the execution of 1D & 2D simulations.
        function initializeShiftingStates(obj)
            obj.pShiftEnabled = false(2,3);
            for i=1:ndims(obj.pArray)
                obj.pShiftEnabled(:, i) = size(obj.pArray, i);
            end
        end
        
%___________________________________________________________________________________________________ initializeBoundingEdges
% Initializes the Edge object that stores and manages the necessary data for operations, such as 
% shifting, on the bounding edges of the computational grid. For this function to work properly, the
% readBoundaryConditions method must be called first and the initial condition for the data array
% must already be set.
        function initializeBoundingEdges(obj)
            obj.edgeStore = Edges(obj.bcModes, obj.pArray, 0.005);
        end
        
%___________________________________________________________________________________________________ readStatics
% Reads the static cell array from the structure provided by as an argument in imogen.m.
%>> statics     Class carrying full information regarding all statics in simulation      class
        function readStatics(obj, statics)
            if isempty(statics); return; end
            
                %--- Flux array case ---%
                if isa(obj,'FluxArray')
                    [SI SV SC] = statics.staticsForVariable(obj.id{1}, obj.component, statics.FLUXL);
                    obj.staticIndices = SI;
                    obj.staticValues  = SV;
                    obj.staticCoeffs  = SC;

                    if isempty(SI); obj.staticActive = false; else; obj.staticActive = true; end
                %--- Primary array case ---%
                else
                    [SI SV SC] = statics.staticsForVariable(obj.id{1}, obj.component, statics.CELLVAR);
                    obj.staticIndices = SI;
                    obj.staticValues  = SV;
                    obj.staticCoeffs  = SC;

                    if isempty(SI); obj.staticActive = false; else; obj.staticActive = true; end
                end

                if (obj.pRunManager.useGPU == true) && obj.staticActive
                    obj.staticLinIndices = GPUdouble(obj.staticIndices(:,1)-1);
                    obj.staticValues = GPUdouble(obj.staticValues);
                    obj.staticCoeffs = GPUdouble(obj.staticCoeffs);
                end

        end
        
    end%PROTECTED
    
%===================================================================================================
    methods (Static = true) %                                                       S T A T I C  [M]
        
%___________________________________________________________________________________________________ zeroIt
% Sets all values below the 1e-12 tolerance to zero on the input array.
% inArray   Input array to be zeroed.                                                   double
        function result = zeroIt(inArray)
           result = inArray .* (inArray > 1e-12);
        end
        
    end%STATIC

end
