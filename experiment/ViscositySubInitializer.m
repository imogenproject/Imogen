classdef ViscositySubInitializer < handle
% Sub initializer class for artificial viscosity variables.
	
%===================================================================================================
	properties (Constant = true, Transient = true) %							C O N S T A N T	 [P]
    end%CONSTANT
	
%===================================================================================================
    properties (SetAccess = public, GetAccess = public) %							P U B L I C  [P]
        type;           % Enumerated artifitical viscosity model type.          string
        linear;         % Linear artificial viscosity coefficient.              double
        quadratic;      % Quadratic artificial viscosity coefficient.           double
    end %PUBLIC

%===================================================================================================
    properties (SetAccess = protected, GetAccess = protected) %				   P R O T E C T E D [P]
    end %PROTECTED
	
	
	
	
	
%===================================================================================================
    methods %																	  G E T / S E T  [M]
	
%___________________________________________________________________________________________________ ViscositySubInitializer
        function obj = ViscositySubInitializer()
            obj.type        = ENUM.ARTIFICIAL_VISCOSITY_NONE;
            obj.linear      = 0;
            obj.quadratic   = 0;
        end
        
    end%GET/SET
	
%===================================================================================================
    methods (Access = public) %														P U B L I C  [M]
	end%PUBLIC
	
%===================================================================================================	
	methods (Access = protected) %											P R O T E C T E D    [M]
	end%PROTECTED
		
%===================================================================================================	
	methods (Static = true) %													  S T A T I C    [M]
	end%PROTECTED
	
end%CLASS
