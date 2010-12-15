classdef JetInitializer < Initializer
% Creates imogen input data arrays for the fluid, magnetic-fluid, and gravity jet tests.
%
% Unique properties for this initializer:
%   direction  % enumerated direction of the jet.                                  str
%   jetMass    % mass value of the jet.                                            double
%   jetMach    % mach speed for the jet.                                           double
%   jetMags    % magnetic field amplitudes for the jet.                            double(3)
%   offset     % index location of the jet on the grid.                            double(3)
%   backMass   % mass value in background fluid.                                   double
%   backMags   % magnetic field amplitudes in background fluid.                    double(3)
        
%===================================================================================================
    properties (Constant = true, Transient = true) %                            C O N S T A N T  [P]
        X = 'x';
        Y = 'y';
        Z = 'z';
        
    end%CONSTANT
    
%===================================================================================================
    properties (SetAccess = public, GetAccess = public) %                           P U B L I C  [P]
        direction;  % enumerated direction of the jet.                                  str
        flip;       % specifies if the jet should be negatively directed.               logical
        jetMass;    % mass value of the jet.                                            double
        jetMach;    % mach speed for the jet.                                           double
        jetMags;    % magnetic field amplitudes for the jet.                            double(3)
        offset;     % index location of the jet on the grid.                            double(3)
        backMass;   % mass value in background fluid.                                   double
        backMags;   % magnetic field amplitudes in background fluid.                    double(3)
    end %PUBLIC

%===================================================================================================
    properties (Dependent = true) %                                            D E P E N D E N T [P]
    end %DEPENDENT
    
%===================================================================================================
    properties (SetAccess = protected, GetAccess = protected) %                P R O T E C T E D [P]
    end %PROTECTED
    
    
    
    
    
%===================================================================================================
    methods %                                                                     G E T / S E T  [M]
        
%___________________________________________________________________________________________________ JetInitializer
        function obj = JetInitializer(input)
            obj                  = obj@Initializer();
            obj.gamma            = 5/3;
            obj.runCode          = 'Jet';
            obj.info             = 'Jet trial.';
            obj.mode.fluid		 = true;
            obj.mode.magnet		 = false;
            obj.mode.gravity	 = false;
            obj.cfl				 = 0.7;
            obj.iterMax          = 250;
            obj.bcMode.x		 = 'trans';
            obj.bcMode.y		 = 'fade';
            obj.bcMode.z         = 'circ';
            obj.activeSlices.xy  = true;
            obj.ppSave.dim2      = 10;
            
            obj.direction       = JetInitializer.X;
            obj.flip            = false;
            obj.jetMass         = 1;
            obj.jetMach         = 1;
            obj.jetMags         = [0 0 0];
            obj.backMass        = 1;
            obj.backMags        = [0 0 0];
            
            obj.operateOnInput(input, [512 256 1]);
        end
        
%___________________________________________________________________________________________________ jetMags
        function result = get.jetMags(obj)
            result = obj.make3D(obj.jetMags, 0);
        end
           
%___________________________________________________________________________________________________ jetMags
        function result = get.backMags(obj)
            result = obj.make3D(obj.backMags, 0);
        end
        
        
    end%GET/SET
    
%===================================================================================================
    methods (Access = public) %                                                     P U B L I C  [M]        
    end%PUBLIC
    
%===================================================================================================    
    methods (Access = protected) %                                          P R O T E C T E D    [M]
        
%___________________________________________________________________________________________________ calculateInitialConditions
        function [mass, mom, ener, mag, statics] = calculateInitialConditions(obj)

            
            if isempty(obj.offset)
                obj.offset = [ceil(obj.grid(1)/10), ceil(obj.grid(2)/2), ceil(obj.grid(3)/2)];
            end
                
            mass    = obj.backMass * ones(obj.grid);
            mom     = zeros([3 obj.grid]);
            mag     = zeros([3 obj.grid]);
            
            %--- Magnetic background ---%
            for i=1:3;    mag(i,:,:,:) = obj.backMags(i)*ones(obj.grid); end
            
            %--- Total energy ---%
            magSquared    = squeeze( sum(mag .* mag, 1) );
            ener          = (mass.^obj.gamma)/(obj.gamma - 1) + 0.5*magSquared;
            
            %--- Static values for the jet ---%
            jetMom            = obj.jetMass*speedFromMach(obj.jetMach, obj.gamma, obj.backMass, ...
                                                            ener(1,1,1), obj.backMags');

            if (obj.flip), jetMom = - jetMom; end

            jetEner = (obj.jetMass^obj.gamma)/(obj.gamma - 1) ...   % internal
                        + 0.5*(jetMom^2)/obj.jetMass ...            % kinetic
                        + 0.5*sum(obj.jetMags .* obj.jetMags, 2);   % magnetic

statics = StaticsInitializer();
             
lind = zeros(obj.grid);
lind(1:end) = 1:prod(obj.grid);

xMin = max(obj.offset(1)-2,1);        xMax = min(obj.offset(1)+2,obj.grid(1));
yMin = max(obj.offset(2)-2,1);        yMax = min(obj.offset(2)+2,obj.grid(2));
zMin = max(obj.offset(3)-2,1);        zMax = min(obj.offset(3)+2,obj.grid(3));

statics.valueSet = {0, obj.jetMass, jetMom, jetEner, obj.jetMags(1), ...
                            obj.jetMags(2), obj.jetMags(3)};

statics.indexSet = {lind(xMin:xMax+1,yMin:yMax,zMin:zMax)};

statics.associateStatics(ENUM.MASS, ENUM.SCALAR, statics.CELLVAR, 1, 2);
%statics.associateStatics(ENUM.MASS, ENUM.SCALAR, statics.FLUXL, 1, 2);
%statics.associateStatics(ENUM.ENER, ENUM.SCALAR, statics.CELLVAR, 1, 4);

%fields = {JetInitializer.X, JetInitializer.Y, JetInitializer.Z};
%for i = 1:3
%	if strcmp(obj.direction, fields{i});    momIndex = 3;
 %       else                                    momIndex = 1;
  %      end
	statics.associateStatics(ENUM.MOM, ENUM.VECTOR(1), statics.CELLVAR, 1, 3);
%end

%statics.associateStatics(ENUM.MAG, ENUM.VECTOR(1), statics.CELLVAR, 1, 5);
%statics.associateStatics(ENUM.MAG, ENUM.VECTOR(2), statics.CELLVAR, 1, 6);
%statics.associateStatics(ENUM.MAG, ENUM.VECTOR(3), statics.CELLVAR, 1, 7);
            
            if obj.mode.magnet;     obj.runCode = [obj.runCode 'Mag'];  end
            if obj.mode.gravity;    obj.runCode = [obj.runCode 'Grav']; end
        end
        
%___________________________________________________________________________________________________ toInfo
        function result = toInfo(obj)
            skips = {'X', 'Y', 'Z'};
            result = toInfo@Initializer(obj, skips);
        end                    
        
    end%PROTECTED
        
%===================================================================================================    
    methods (Static = true) %                                                     S T A T I C    [M]
    end
end%CLASS
