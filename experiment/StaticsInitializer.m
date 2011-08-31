classdef StaticsInitializer < handle
%___________________________________________________________________________________________________ 
%===================================================================================================
    properties (Constant = true, Transient = true) %                     C O N S T A N T         [P]
        CELLVAR = 0;
        FLUXL   = 1;
        FLUXR   = 2;
        FLUXALL = 3;
    end%CONSTANT
        
%===================================================================================================
    properties (SetAccess = public, GetAccess = public) %                           P U B L I C  [P]

        indexSet; % Arrays in indices to set statically            Cell[N]
        valueSet; % Arrays of values to set statically             Cell[N]
        coeffSet; % Arrays of coefficients to set [must match dimensions of corresponding value set]

        arrayStatics; % Cell array with one cell per simulation var Cell[8];
        % WARNING: THIS MUST BE THE SAME SIZE AS THE NUMBER OF SIMULATION VARIABLES

    end %PUBLIC

%===================================================================================================
    properties (SetAccess = protected, GetAccess = protected) %                P R O T E C T E D [P]

        readyForReadout; % Set if prepareStaticsForSimulation() has been called and
                         % no new associations have been created since

    end %PROTECTED
        
%===================================================================================================
    methods %                                                                     G E T / S E T  [M]

    end%GET/SET
        
%===================================================================================================
    methods (Access = public) %                                                     P U B L I C  [M]

        function obj = StaticsInitializer()
            obj.arrayStatics = cell(8,1); % Create one arrayStatics for every variable
            for x = 1:8
                obj.arrayStatics{x} = struct('arrayField',[], 'indexId',[], 'valueId',[], 'coeffId',[]);
                % Create structs to associate pairs of index sets and values with the primary & flux
                % arrays of each simulation variable
            end
        end

        function [indices values coeffs] = staticsForVariable(obj, varId, component, fieldId)
            indices = [];
            values = [];
            coeffs = [];

            if ~obj.readyForReadout; obj.prepareStaticsForSimulation(); end

            partmap = obj.mapVaridToIdx(varId, component);
            if partmap == 0; return; end

            AS = obj.arrayStatics{partmap};

            for x = 1:numel(AS.arrayField) % For each static defined for this variable
                if AS.arrayField(x) == fieldId % If it applies to the requested field

                    newIdx   = obj.indexSet{AS.indexId(x)};
                    newVal   = obj.valueSet{AS.valueId(x)};

                    if AS.coeffId(x) == 0
                        newCoeff = 1;
                    else
                        newCoeff = obj.coeffSet{AS.coeffId(x)};
                    end

                    % Expand array-scalar pairs to array-array pairs; This can be done
                    if numel(newVal) == 1; newVal   = newVal   * ones(size(newIdx,1),1); end
                    if numel(newCoeff)==1; newCoeff = newCoeff * ones(size(newIdx,1),1); end

                    % Fail if nonequally sized arrays are paired; This cannot be done
                    if size(newVal,1) ~= size(newIdx,1)
                    error(sprintf('Unrecoverable error preparing statics; numel(index set %i) = %i but numel(value set %i) = %i.\n', x, size(newIdx,1), x, numel(newVal)));
                    end

                    indices = [indices; newIdx]; % cat indices
                    values  = [values ; newVal]; % cat values
                    coeffs  = [coeffs ; newCoeff]; % cat fade coefficients
                end
            end

        end

        % This function prepares statics for injection into array statics by reshaping for concatenation
        function prepareStaticsForSimulation(obj)
            % Reshape them to be Nx1 arrays so we can cat using [u; v]
            for x = 1:numel(obj.indexSet)
                % Reshape them to be Nx1
%                obj.indexSet{x} = reshape(obj.indexSet{x}, [numel(obj.indexSet{x}) 1]);
%                obj.valueSet{x} = reshape(obj.valueSet{x}, [numel(obj.valueSet{x}) 1]);
            end

            obj.readyForReadout = 1;

        end

        % Adds a pair of statics
        function addStatics(obj, indices, values, coeffs)
            obj.indexSet{end+1} = indices;
            obj.valueSet{end+1} = values;
            if nargin == 3
                obj.coeffSet{end+1} = 1;
            else
                obj.coeffSet{end+1} = coeffs;
            end
        end

        % Maps a set of indices, values and fade rate coefficients to a variable
        function associateStatics(obj, varID, component, fieldID, indexNum, valueNum, coeffNum)
            vmap = obj.mapVaridToIdx(varID, component);

            obj.arrayStatics{vmap}.arrayField(end+1) = fieldID;
            obj.arrayStatics{vmap}.indexId(end+1)    = indexNum;
            obj.arrayStatics{vmap}.valueId(end+1)    = valueNum;
            if nargin == 6
                obj.arrayStatics{vmap}.coeffId(end+1) = 0;
            else
                obj.arrayStatics{vmap}.coeffId(end+1)    = coeffNum;
            end

            obj.readyForReadout = 0;
        end

        %%% === Utility functions === %%%

        function setFluid_allConstantBC(obj, mass, ener, mom, facenumber)
            obj.setConstantBC(ENUM.MASS, ENUM.SCALAR,   obj.CELLVAR, mass, facenumber);
            obj.setConstantBC(ENUM.ENER, ENUM.SCALAR,   obj.CELLVAR, ener, facenumber);
            obj.setConstantBC(ENUM.MOM, ENUM.VECTOR(1), obj.CELLVAR, squeeze(mom(1,:,:,:)), facenumber);
            obj.setConstantBC(ENUM.MOM, ENUM.VECTOR(2), obj.CELLVAR, squeeze(mom(2,:,:,:)), facenumber);
            if size(mom,4) > 1
                obj.setConstantBC(ENUM.MOM, ENUM.VECTOR(3), ENUM.CELLVAR, squeeze(mom(3,:,:,:)), facenumber);
            end

        end

        function setMag_allConstantBC(obj, mag, facenumber)
            obj.setConstantBC(ENUM.MAG, ENUM.VECTOR(1), obj.CELLVAR, squeeze(mag(1,:,:,:)), facenumber);
            obj.setConstantBC(ENUM.MAG, ENUM.VECTOR(2), obj.CELLVAR, squeeze(mag(2,:,:,:)), facenumber);
            if size(mag,4) > 1
                obj.setConstantBC(ENUM.MAG, ENUM.VECTOR(3), ENUM.CELLVAR, squeeze(mag(3,:,:,:)), facenumber);
            end
        end

        function setFluid_allFadeBC(obj, mass, ener, mom, facenumber, bcinf)
            obj.setFadeBC(ENUM.MASS, ENUM.SCALAR,   obj.CELLVAR, mass, facenumber, bcinf);
            obj.setFadeBC(ENUM.ENER, ENUM.SCALAR,   obj.CELLVAR, ener, facenumber, bcinf);
            obj.setFadeBC(ENUM.MOM, ENUM.VECTOR(1), obj.CELLVAR, squeeze(mom(1,:,:,:)), facenumber, bcinf);
            obj.setFadeBC(ENUM.MOM, ENUM.VECTOR(2), obj.CELLVAR, squeeze(mom(2,:,:,:)), facenumber, bcinf);
            if size(mom,4) > 1
                obj.setFadeBC(ENUM.MOM, ENUM.VECTOR(3), ENUM.CELLVAR, squeeze(mom(3,:,:,:)), facenumber, bcinf);
            end
        end

        function setMag_allFadeBC(obj, mag, facenumber, bcinf)
            obj.setFadeBC(ENUM.MAG, ENUM.VECTOR(1), obj.CELLVAR, squeeze(mag(1,:,:,:)), facenumber, bcinf);
            obj.setFadeBC(ENUM.MAG, ENUM.VECTOR(2), obj.CELLVAR, squeeze(mag(2,:,:,:)), facenumber, bcinf);
            if size(mag,4) > 1
                obj.setFadeBC(ENUM.MAG, ENUM.VECTOR(3), ENUM.CELLVAR, squeeze(mag(3,:,:,:)), facenumber, bcinf);
            end
        end

        function setConstantBC(obj, varID, component, fieldID, array, facenumber)
            vmap = obj.mapVaridToIdx(varID, component);

            xset=[]; yset=[]; zset=[];

            switch facenumber
                case 1; xset=1:2;
                        yset=1:size(array,2); zset=1:size(array,3); % minus X
                case 2; xset=(size(array,1)-1):size(array,1);
                        yset=1:size(array,2); zset=1:size(array,3);% plux  X

                case 3; xset=1:size(array,1); % minus Y
                        yset=1:2; zset=1:size(array,3);
                case 4; xset=1:size(array,1); % plus  Y
                        yset=(size(array,2)-1):size(array,2); zset=1:size(array,3);

                case 5; xset=1:size(array,1); yset=1:size(array,2); % minus Z
                        zset=1:2;
                case 6; xset=1:size(array,1); yset=1:size(array,2); % plus  Z
                        zset=(size(array,3)-1):size(array,3);
            end

            inds = indexSet(size(array), xset, yset, zset);

            obj.addStatics(inds, array(inds(:,1)));

            obj.associateStatics(varID, component, fieldID, numel(obj.indexSet), numel(obj.valueSet), numel(obj.coeffSet));
        end

        % DO NOT USE YET
        function setFadeBC(obj, varID, component, fieldID, array, facenumber, bcInfinity)
            vmap = obj.mapVaridToIdx(varID, component);

            xset=[]; yset=[]; zset=[];
            coeff = [];
            AS = size(array); if numel(AS) == 2; AS(3) = 1; end
            switch facenumber
                case 1; xset=1:bcInfinity;
                        yset=1:AS(2); zset=1:AS(3); % minus X
                        coeff = ndgrid(1:bcInfinity, yset, zset);
                case 2; xset=AS(1):-1:(AS(1)-bcInfinity+1);
                        yset=1:AS(2); zset=1:AS(3);% plux  X
                        coeff = ndgrid(bcInfinity:-1:1, yset, zset);

                case 3; xset=1:AS(1); % minus Y
                        yset=1:bcInfinity; zset=1:AS(3);
                        [drop1 coeff drop2] = ndgrid(xset, 1:bcInfinity, zset);
                case 4; xset=1:AS(1); % plus  Y
                        yset=AS(2):-1:(AS(2)-bcInfinity+1); zset=1:AS(3);
                        [drop1 coeff drop2] = ndgrid(xset, bcInfinity:-1:1, zset);

                case 5; xset=1:AS(1); yset=1:AS(2); % minus Z
                        zset=1:bcInfinity;
                        [drop1 drop2 coeff] = ndgrid(xset, yset, 1:bcInfinity);
                case 6; xset=1:AS(1); yset=1:AS(2); % plus  Z
                        zset=AS(3):-1:(AS(3)-bcInfinity+1);
                        [drop1 drop2 coeff] = ndgrid(xset, yset, bcInfinity:-1:1);
            end

            inds = indexSet(size(array), xset, yset, zset);
            coeff= pchip([0 ceil(bcInfinity/4) round(bcInfinity/2) (bcInfinity-1) bcInfinity], [1 1 .02 0 0], coeff); 
            
            obj.addStatics(inds, array(inds(:,1)), coeff(:));

            obj.associateStatics(varID, component, fieldID, numel(obj.indexSet), numel(obj.valueSet), numel(obj.coeffSet));
        end

    end%PUBLIC
        
%===================================================================================================        
    methods (Access = protected) %                                          P R O T E C T E D    [M]

    end%PROTECTED
                
%===================================================================================================        
    methods (Static = true) %                                                     S T A T I C    [M]

        function result = mapVaridToIdx(varId, component)
            if strcmp(varId,ENUM.MASS); result = 1; return; end
            if strcmp(varId,ENUM.ENER); result = 2; return; end
            if strcmp(varId,ENUM.MOM); result = 2+component; return; end % 3-4-5
            if strcmp(varId,ENUM.MAG); result = 5+component; return; end % 6-7-8

            result = 0;
            return;
        end

    end%PROTECTED
        
end%CLASS
