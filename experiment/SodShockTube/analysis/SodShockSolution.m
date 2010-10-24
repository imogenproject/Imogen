classdef SodShockSolution < handle
% Class annotation template for creating new classes.
	
	
%___________________________________________________________________________________________________ 
	
%===================================================================================================
	properties (Constant = true, Transient = true) %							C O N S T A N T	 [P]
        GAMMA       = 1.4;
    end%CONSTANT
	
%===================================================================================================
    properties (SetAccess = public, GetAccess = public) %							P U B L I C  [P]
        time;        
        x;
        mass;
        soundSpeed;
        pressure;
        velocity;
        energy;
    end %PUBLIC

%===================================================================================================
    properties (SetAccess = protected, GetAccess = protected) %				   P R O T E C T E D [P]
        pTubeLength;        % Length of the tube.                                       double
    end %PROTECTED
	
	
	
%===================================================================================================
    methods %																	  G E T / S E T  [M]

%___________________________________________________________________________________________________ SodShockSolution
        function obj = SodShockSolution(resolution, time)
            
            obj.time                    = time;
            g                           = obj.GAMMA;
            x0                          = floor(resolution/2);
            
            obj.mass                    = zeros(1, resolution);
            obj.mass(1:x0)              = 1;
            obj.mass(x0:end)            = 0.125;
            
            obj.pressure                = zeros(1, resolution);
            obj.pressure(1:x0)          = 1;
            obj.pressure(x0:end)        = 0.1;

            obj.velocity                = zeros(1, resolution);
            obj.soundSpeed              = sqrt(g*obj.pressure./obj.mass);
            
            %--- Post Shock ---%
            postPressure                = obj.findPostShockPressure();
            postVelocity                = obj.calculatePostShockVelocity(postPressure);
            postMass                    = obj.calculatePostShockMass(postPressure);
            
            shockSpeed                  = obj.calculateShockSpeed(postMass, postVelocity);
            postContactMass             = obj.calculatePostContactMass(postPressure);
                        
            %--- Rarefaction ---%
            rarefactionSoundSpeed       = obj.calculateRarefactionSoundSpeed(0.5, time);
            rarefactionVelocity         = obj.calculateRarefactionVelocity(0.5, time);
            rarefactionMass             = obj.calculateRarefactionMass(rarefactionSoundSpeed);
            rarefactionPressure         = obj.calculateRarefactionPressure(rarefactionMass);
            
            %--- Find Positions ---%
            x1                          = ceil((0.5 - obj.soundSpeed(1)*time)*resolution);
            x2                          = find(rarefactionMass <= postContactMass, true, 'first');
            x3                          = floor((0.5 + postVelocity*time)*resolution);
            x4                          = floor((0.5 + shockSpeed*time)*resolution);
            
            obj.mass(x1:x2)             = rarefactionMass(x1:x2);
            obj.mass(x2:x3)             = postContactMass;
            obj.mass(x3:x4)             = postMass;
            
            obj.pressure(x1:x2)         = rarefactionPressure(x1:x2);
            obj.pressure(x2:x4)         = postPressure;
            
            obj.velocity(x1:x2)         = rarefactionVelocity(x1:x2);
            obj.velocity(x2:x4)         = postVelocity;
            
            obj.soundSpeed              = sqrt(g*obj.pressure./obj.mass);
            obj.energy                  = obj.pressure./(g - 1) ...
                                            + 0.5*obj.mass.*obj.velocity.^2;
        end
        
	end%GET/SET
	
%===================================================================================================
    methods (Access = public) %														P U B L I C  [M]
	end%PUBLIC
	
%===================================================================================================	
	methods (Access = protected) %											P R O T E C T E D    [M]
        
%___________________________________________________________________________________________________ findPostShockPressure
        function result = findPostShockPressure(obj)
            g = obj.GAMMA;
            G = (g-1)/(g+1);
            b = (g-1)/(2*g);
            
            m1 = obj.mass(1);
            m2 = obj.mass(end);
            p1 = obj.pressure(1);
            p2 = obj.pressure(end);
            
            func = @(p)...
                    (p1^b - p^b)*sqrt((1-G^2)*p1^(1/g)/(G^2*m1)) ...
                    - (p - p2)*sqrt((1-G)/(m2*(p+G*p2)));
            result = fzero(func, 0.3);
        end
        
%___________________________________________________________________________________________________ calculatePostShockVelocity
        function result = calculatePostShockVelocity(obj, postPressure)
            g      = obj.GAMMA;
            result = 2*sqrt(g)/(g - 1)*(1 - postPressure^((g - 1)/(2*g)));
        end
        
%___________________________________________________________________________________________________ calculatePostShockMass
        function result = calculatePostShockMass(obj, postPressure)
            g      = obj.GAMMA;
            G      = (g-1)/(g+1);
            m2     = obj.mass(end);
            p2     = obj.pressure(end);
            
            result = m2*((postPressure/p2) + G) / (1 + G*(postPressure/p2));
        end
        
%___________________________________________________________________________________________________ calculateShockSpeed
        function result = calculateShockSpeed(obj, postMass, postVelocity)
            m2     = obj.mass(end);
            
            result = postVelocity*(postMass/m2) / ((postMass/m2) - 1);  
        end
        
%___________________________________________________________________________________________________ calculatePostContactMass
        function result = calculatePostContactMass(obj, postPressure)
            g      = obj.GAMMA;
            p1     = obj.pressure(1);
            
            result = obj.mass(1)*(postPressure/p1)^(1/g);
        end
        
%___________________________________________________________________________________________________ calculateRarefactionSoundSpeed
        function result = calculateRarefactionSoundSpeed(obj, x0, time)
            g         = obj.GAMMA;
            G         = (g-1)/(g+1);
            c1        = obj.soundSpeed(1);
            positions = linspace(0, 1, length(obj.mass));
            
            result    = G*((x0-positions)./time) + (1 - G)*c1;
        end
        
%___________________________________________________________________________________________________ calculateRarefactionVelocity
        function result = calculateRarefactionVelocity(obj, x0, time)
            g         = obj.GAMMA;
            G         = (g-1)/(g+1);
            c1        = obj.soundSpeed(1);
            positions = linspace(0, 1, length(obj.soundSpeed));
            
            result    = (1 - G)*((positions-x0)/time + c1); 
        end
        
%___________________________________________________________________________________________________ calculateRarefactionMass
        function result = calculateRarefactionMass(obj, rarefactionSoundSpeed)
            g      = obj.GAMMA;
            m1     = obj.mass(1);
            c1     = obj.soundSpeed(1);
            
            result = m1.*(rarefactionSoundSpeed./c1).^(2/(g-1)); 
        end
        
%___________________________________________________________________________________________________ calculateRarefactionPressure
        function result = calculateRarefactionPressure(obj, rarefactionMass)
            g      = obj.GAMMA;
            m1     = obj.mass(1);
            p1     = obj.pressure(1);
            
            result = p1*(rarefactionMass./m1).^g;
        end
        
	end%PROTECTED
		
%===================================================================================================	
	methods (Static = true) %													  S T A T I C    [M]
	end%PROTECTED
	
end%CLASS