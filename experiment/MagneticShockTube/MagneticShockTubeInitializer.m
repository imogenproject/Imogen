classdef MagneticShockTubeInitializer < Initializer
% Creates initial conditions for a magnetic shock tube simulation. This is a canonical test for 
% magnetichydrodynamic evolution in MHD codes. While there is no analytical solution as is the case
% for the Sod shock tube problem, the evolution of the magnetic shock tube has been widely analyzed
% by many different MHD codes, resulting in a consensus on the long-term evolutionary behavior.
%
% Unique properties for this initializer:
%   direction;      % Enumerated spatial orientation of the shock wave                     str
    
        
%===================================================================================================
    properties (Constant = true, Transient = true) %                            C O N S T A N T  [P]
        X = 'X';
        Y = 'Y';
        Z = 'Z';
    end%CONSTANT
    
%===================================================================================================
    properties (SetAccess = public, GetAccess = public) %                           P U B L I C  [P]
        direction;      % Enumerated spatial orientation of the shock wave.         str
        xField;         % Specifies if the x-magnetic field is non-zero.            logical
        yField;         % Specifies if the y-magnetic field is non-zero.            logical
        zField;         % Specifies if the z-magnetic field is non-zero.            logical
    end %PUBLIC

%===================================================================================================
    properties (Dependent = true) %                                            D E P E N D E N T [P]
    end %DEPENDENT
    
%===================================================================================================
    properties (SetAccess = protected, GetAccess = protected) %                P R O T E C T E D [P]
    end %PROTECTED
    
    
    
    
    
%===================================================================================================
    methods %                                                                     G E T / S E T  [M]
        
%___________________________________________________________________________________________________ MagneticShockTubeInitializer
        function obj = MagneticShockTubeInitializer(input)
            obj                  = obj@Initializer();
            obj.gamma            = 1.4;
            obj.runCode          = 'MagST';
            obj.info             = 'Magnetic Shock tube trial.';
            obj.mode.fluid		 = true;
            obj.mode.magnet		 = true;
            obj.mode.gravity	 = false;
            obj.cfl				 = 0.35;
            obj.iterMax          = 250;
            obj.bcMode.x		 = 'fade';
            obj.bcMode.y		 = 'circ';
            obj.bcMode.z         = 'circ';
            obj.activeSlices.x   = true;
            obj.activeSlices.xyz = true;
            obj.ppSave.dim1      = 10;
            obj.ppSave.dim3      = 25;
            
            obj.direction       = 'X';
            obj.xField          = true;
            obj.yField          = true;
            obj.zField          = false;
            
            obj.operateOnInput(input, [1024 3 3]);
        end
               
        
    end%GET/SET
    
%===================================================================================================
    methods (Access = public) %                                                     P U B L I C  [M]
    end%PUBLIC
    
%===================================================================================================    
    methods (Access = protected) %                                          P R O T E C T E D    [M]
        
%___________________________________________________________________________________________________ calculateInitialConditions
        function [mass, mom, ener, mag, statics, run] = calculateInitialConditions(obj)
        
            %--- Initialization ---%
            statics = [];
            half    = round(obj.grid/2);
            indeces = cell(1,3);
            for i=1:3;  indeces{i} = 1:obj.grid(i); end

            %--- Set values for output arrays ---%
            mass    = ones(obj.grid);
            mom     = zeros([3 obj.grid]);
            mag     = 5/sqrt(4*pi)*ones([3 obj.grid]); %Start with non-zero mags in all directions

            %--- Zero magnetic field components according to mags argument values ---%
            if ~obj.xField;     mag(1,:,:,:) = zeros([1 obj.grid]); end
            if ~obj.yField;     mag(2,:,:,:) = zeros([1 obj.grid]); end
            if ~obj.zField;     mag(3,:,:,:) = zeros([1 obj.grid]); end

            %--- Adjust momentum values for each side of the shock ---%
            direct = {'X','Y','Z'};
            for i=1:3
                if strcmpi(obj.direction, direct{i})
                    indeces{i} = half(i):obj.grid(i);
                    mom(i,:,:,:) = 10 * ones([1 obj.grid]);
                    mom(i,indeces{:}) = -1.0 * mom(1,indeces{:});
                    obj.runCode = [obj.runCode direct{i}];
                    break;
                end
            end
    
            %--- Calculate energy ---%
            %       This problem assumes the LHS of the shock has a pressure of 1 and the RHS has
            %       a pressure of 20. Internal energies are calculated based on these values.
            ener             = 0.5*squeeze(sum(mom.*mom,1))./mass ...   % kinetic energy
                             + 0.5*squeeze(sum(mag.*mag,1));            % magnetic energy
            ener(indeces{:}) = ener(indeces{:}) + 1/(obj.gamma - 1)*ones(size(ener(indeces{:})));
            indeces{i}       = 1:(indeces{i}(1)-1);
            ener(indeces{:}) = ener(indeces{:}) + 20/(obj.gamma - 1)*ones(size(ener(indeces{:})));
        
            run = obj.getRunSettings();
        end
        
    end%PROTECTED
        
%===================================================================================================    
    methods (Static = true) %                                                     S T A T I C    [M]
    end
end%CLASS
