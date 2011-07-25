classdef SodShockTubeInitializer < Initializer
% Creates initial conditions for a Sod shock tube simulation. This is a canonical test of the 
% hydrodynamic portion of MHD codes as the solution can be determined analytically for comparison
% to metric the functionality of a code. Imogen has been tested thoroughly with this simulation and,
% as such, it is an excellent tool to verify the continued operation of the fluid routines.
%
% Unique properties for this initializer:
%   direction       % Enumerated spatial orientation of the shock wave                      str
%   shockAngle      % Off-axis angle for the shockwave in degrees.                          double
    
        
%===================================================================================================
    properties (Constant = true, Transient = true) %                            C O N S T A N T  [P]
        X = 'X';
        Y = 'Y';
        Z = 'Z';
    end%CONSTANT
    
%===================================================================================================
    properties (SetAccess = public, GetAccess = public) %                           P U B L I C  [P]
        direction;      % Enumerated spatial orientation of the shock wave       	str
        shockAngle;     % Angle of the shock axis                                   double
    end %PUBLIC

%===================================================================================================
    properties (Dependent = true) %                                            D E P E N D E N T [P]
    end %DEPENDENT
    
%===================================================================================================
    properties (SetAccess = protected, GetAccess = protected) %                P R O T E C T E D [P]
    end %PROTECTED
    
    
    
    
    
%===================================================================================================
    methods %                                                                     G E T / S E T  [M]
        
%___________________________________________________________________________________________________ SodShockTubeInitializer
        function obj = SodShockTubeInitializer(input)
            obj                  = obj@Initializer();
            obj.gamma            = 1.4;
            obj.runCode          = 'SodST';
            obj.info             = 'Sod shock tube trial.';
            obj.mode.fluid		 = true;
            obj.mode.magnet		 = false;
            obj.mode.gravity	 = false;
            obj.cfl				 = 0.7;
            obj.iterMax          = 150;
            obj.ppSave.dim1      = 10;
            obj.ppSave.dim3      = 25;
            
            obj.shockAngle       = 0;
            obj.direction        = 'X';
            
            obj.bcMode.x         = ENUM.BCMODE_CONST;
            obj.bcMode.y         = ENUM.BCMODE_CONST;
            obj.bcMode.z         = ENUM.BCMODE_CONST;
            
            obj.operateOnInput(input, [1024, 4, 4]);
        end
        
        
    end%GET/SET
    
%===================================================================================================
    methods (Access = public) %                                                     P U B L I C  [M]        
    end%PUBLIC
    
%===================================================================================================    
    methods (Access = protected) %                                          P R O T E C T E D    [M]
        
%___________________________________________________________________________________________________ calculateInitialConditions
        function [mass, mom, ener, mag, statics] = calculateInitialConditions(obj)
        
            %--- Initialization ---%
            obj.runCode           = [obj.runCode upper(obj.direction)];
            statics               = []; % No statics used in this problem
            half                  = round(obj.grid/2);
            indices               = cell(1,3);
            for i=1:3
                indices{i} = 1:obj.grid(i);  
            end

            %--- Set array values ---%
            mass                  = ones(obj.grid);
            mom                   = zeros([3, obj.grid]);
            mag                   = zeros([3, obj.grid]);
            ener                  = ones(obj.grid) / (obj.gamma - 1);

            %--- Set shock array values according to flux direction ---%
            direct                = {'x', 'y', 'z'};
            i                     = find(strcmpi(obj.direction, direct), true);
            j                     = mod(i,2) + 1;
            
            for n=1:obj.grid(j)
                adjacentLen       = -half(j) + n - 1;
                lowerBound        = floor(half(i) ...
                                    - adjacentLen*tand(obj.shockAngle));
                lowerBound        = min(obj.grid(i), max(1, lowerBound));

                indices{i}        = lowerBound:obj.grid(i);
                indices{j}        = n;

                mass(indices{:})  = 0.125*mass(indices{:});
                ener(indices{:})  = 0.1*ener(indices{:});
            end
            
            %--- Adjust Cell Spacing ---%
            %       Problem is normalized so that the length from one end to the other of the shock
            %       tube is 1 unit length, no matter what the resolution. If the shock tube is 
            %       angled, then the normalization is scaled so that the end to end length of the 
            %       shock tube is unit length when cut down the center of tube along the shock
            %       normal.
            obj.dGrid             = 1./obj.grid;
            
            if obj.shockAngle > 0
                angle             = obj.shockAngle;
                criticalAngle     = atan(obj.grid(j)/obj.grid(i));
                if angle <= criticalAngle 
                    scale         = 1/(obj.grid(i)/cosd(angle));
                else
                    scale         = 1/(obj.grid(j)/sind(angle));
                end
                obj.dGrid(i)      = scale;
                obj.dGrid(j)      = scale;
            end
            
            %--- Determine the default slices to save ---%
            if ~obj.saveSlicesSpecified
                obj.activeSlices.(direct{i}) = true;
                
                if obj.shockAngle > 0
                    switch i+j
                        case 3     
                            obj.activeSlices.xy = true;
                        case 5
                            obj.activeSlices.yz = true;
                    end
                end
                
                obj.activeSlices.xyz = true;
            end
            
        end
        
    end%PROTECTED
        
%===================================================================================================    
    methods (Static = true) %                                                     S T A T I C    [M]
    end
end%CLASS
