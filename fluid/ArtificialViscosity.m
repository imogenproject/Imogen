classdef ArtificialViscosity < handle
% Encapsulation class holding all of the artificial viscosity for a simulation.
    
%===================================================================================================
    properties (Constant = true, Transient = true) %                            C O N S T A N T  [P]
    end%CONSTANT
    
%===================================================================================================
    properties (SetAccess = public, GetAccess = public, Transient = true) %         P U B L I C  [P]
        type;                     % Enumerated artificial viscosity type.               string
        linearViscousStrength;    % Coefficient for linear artificial viscosity.        double
        quadraticViscousStrength; % Coefficient for quadratic artificial viscosity.     double
        solve;                    % Handle to viscosity function for run.               @func
    end%PUBLIC
    
%===================================================================================================
    properties (SetAccess = public, GetAccess = private) %                        P R I V A T E  [P]
        pRho0;                    % Storage for density dependence multiplier.          double
    end %PRIVATE
    
    
    
    
    
    
%===================================================================================================
    methods %                                                                     G E T / S E T  [M]
        
%___________________________________________________________________________________________________ ArtificialViscosity
% Creates a new FluidManager instance.
        function obj = ArtificialViscosity() 
            obj.type                     = ENUM.ARTIFICIAL_VISCOSITY_NONE;
            obj.linearViscousStrength    = 0;
            obj.quadraticViscousStrength = 0;
            obj.pRho0                    = 0;
        end
        
    end%GET/SET
    
    
%===================================================================================================
    methods (Access = public) %                                                     P U B L I C  [M]

%___________________________________________________________________________________________________ preliminary
        function preliminary(obj)
            
            switch (obj.type)
                case ENUM.ARTIFICIAL_VISCOSITY_NONE
                    obj.solve = @obj.emptyViscosity;
                case ENUM.ARTIFICIAL_VISCOSITY_NEUMANN_RICHTMYER
                    obj.solve = @obj.neumannRichtmyerViscosity;
                case ENUM.ARTIFICIAL_VISCOSITY_CARAMANA_SHASHKOV_WHALEN
                	obj.solve = @obj.caramanaShashkovWhalenViscosity;
                otherwise
                    obj.solve = @obj.emptyViscosity;
            end
        end
        
    end%PUBLIC
    
    
%===================================================================================================    
    methods (Access = private) %                                                P R I V A T E    [M]

%___________________________________________________________________________________________________  emptyViscosity
        function result = emptyViscosity(obj, run, mass, mom, ener, mag, soundSpeed, direction)
            result = 0;
        end
        
%___________________________________________________________________________________________________ neumannRichtmyerViscosity
% Calculates the Neumann & Richtmyer artificial viscosity term, which smooths out shocks and other
% high frequency effects as a quasi-non-local addition to the pressure.
        function result = neumannRichtmyerViscosity(obj,run,mass,mom,ener,mag,soundSpeed,direction)
            viscousMagnitude = min((mom(direction).shift(direction, 1) - mom(direction).array), 0);
            result = viscousMagnitude .* (obj.linearViscousStrength .* soundSpeed ...
                     + obj.quadraticViscousStrength .* viscousMagnitude ./ mass.array);
        end
        
%___________________________________________________________________________________________________ caramanaShashkovWhalenViscosity
% Calculates artificial viscosity term based on the 1998 paper by Caramana, Shashkov, and Whalen.
		function result = caramanaShashkovWhalenViscosity(obj,run,mass,mom,ener,mag,soundSpeed,direction)				

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Uncomment this section to include averaging.
	%
			% Find minimum of left and right sound speeds.
			massB.array = mass.shift(direction,-1);
			massC.array = mass.shift(direction,1);
%{
			enerB.array = ener.shift(direction,-1);
			enerC.array = ener.shift(direction,1);
			for i=1:3
				momB(i).array = mom(i).shift(direction,-1);
				momC(i).array = mom(i).shift(direction,1);
				magB(i).cellMag.array = mag(i).cellMag.shift(direction,-1);
				magC(i).cellMag.array = mag(i).cellMag.shift(direction,1);
			end
			[soundSpdB, aux] = pressure('speed', run, massB, momB, enerB, magB);
			[soundSpdC, aux] = pressure('speed', run, massC, momC, enerC, magC);
			soundSpeed = min(soundSpdB, soundSpdC);

			% Get velocities of required zones.
			v0 = mom(direction).shift(direction,-2) ./ mass.shift(direction,-2);			
			v1 = momB(direction).array ./ massB.array;
			v2 = momC(direction).array ./ massC.array;
			v3 = mom(direction).shift(direction,2) ./ mass.shift(direction,2);
%}
			% Find average of left and right densities.
			rho  = 2 .* massC.array .* massB.array ./ (massC.array + massB.array);
			soundSpeed = ArtificialViscosity.cswSoundSpeed(run, mass, mom, ener, mag, direction);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Uncomment this section if not including averaging.
	%
 			rho = mass.array;
%}
 			% Get velocities of required zones.
			v0 = mom(direction).shift(direction,-2) ./ mass.shift(direction,-2);			
			v1 = mom(direction).shift(direction,-1) ./ mass.shift(direction,-1);
			v2 = mom(direction).shift(direction,1) ./ mass.shift(direction,1);
			v3 = mom(direction).shift(direction,2) ./ mass.shift(direction,2);	
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			% Calculate limiter using velocity derivative ratios.
			dv = (v2 - v1);
			rL = (v3 - v2) ./ (2 * dv);
			rR = (v1 - v0) ./ (2 * dv);			
			limiter = max(0,min(min(.5 .* (rL + rR), min(2 .* rL, 2 .* rR)), 1));
			
			% Assign viscosity coefficients.
			c1 = obj.linearViscousStrength;
			c2 = obj.quadraticViscousStrength;
%{
            % Make coefficients dependent on mass density.
            if( obj.pRho0 == 0 )
                obj.pRho0 = max( mass.array(:) ) / 1000; % 1/1000 max mass at start
            end
            c1 = 1 + c1 * exp( -( mass.array / obj.pRho0 )); % also try using avg. mass at current
            c2 = 1 + c2 * exp( -( mass.array / obj.pRho0 ));
%}
			% Apply viscosity where velocity derivative is negative.
			result = (dv < 0) .* rho .* (c2 .* (run.GAMMA + 1) .* abs(dv) ./ 4 +...
					 sqrt(c2 ^ 2 .* ((run.GAMMA + 1) / 4) ^ 2 .* dv .^ 2 + c1 ^ 2 .*...
					 soundSpeed .^ 2)) .* (1 - limiter) .* abs(dv) ;
    	end
        
    end%PROTECTED
        
%===================================================================================================    
    methods (Static = true) %                                                     S T A T I C    [M]

%___________________________________________________________________________________________________
% NEED TO REPLACE WITH PRESSURE() FUNCTION CALL, BUT SO FAR IT FAILS.
		function result = cswSoundSpeed(run, mass, mom, ener, mag, direction)
			GAMMA = run.GAMMA;
			
			% Prepare the velocity squared array
			velSquared = mom(1).shift(direction,1) .^ 2 + mom(2).shift(direction,1) .^ 2 ...
						  + mom(3).shift(direction,1) .^ 2 ./ mass.shift(direction,1) .^ 2;
			
			% Prepare the cell-centered magnet squared array
			magSquared = mag(1).cellMag.shift(direction,1) .^2 + mag(2).cellMag.shift(direction,1)...
						 .^2 + mag(3).cellMag.shift(direction,1) .^2;
			
			%Calculate the fluid pressure
			sndSpdC = (GAMMA - 1.0)*(ener.shift(direction,1) - 0.5*mass.shift(direction,1) .* velSquared);
			
			% Calculate the sound speed
			sndSpdC = sndSpdC - (GAMMA - 1.0)*0.5*magSquared;
			sndSpdC = sqrt(abs( (GAMMA*sndSpdC + 2.0*magSquared) ./ mass.shift(direction,1) ));
			
			% Prepare the velocity squared array
			velSquared = mom(1).shift(direction,-1) .^ 2 + mom(2).shift(direction,-1) .^ 2 ...
						  + mom(3).shift(direction,-1) .^ 2 ./ mass.shift(direction,-1) .^ 2;
			
			% Prepare the cell-centered magnet squared array
			magSquared = mag(1).cellMag.shift(direction,-1) .^2 + mag(2).cellMag.shift(direction,-1)...
						 .^2 + mag(3).cellMag.shift(direction,-1) .^2;
			
			%Calculate the fluid pressure
			sndSpdB = (GAMMA - 1.0)*(ener.shift(direction,-1) - 0.5*mass.shift(direction,-1) .* velSquared);
			
			% Calculate the sound speed
			sndSpdB = sndSpdB - (GAMMA - 1.0)*0.5*magSquared;
			sndSpdB = sqrt(abs( (GAMMA*sndSpdB + 2.0*magSquared) ./ mass.shift(direction,-1) ));
			result = min(sndSpdB, sndSpdC);
			result = max( result, 0);
		end

    end%STATIC
    
end%CLASS
        
