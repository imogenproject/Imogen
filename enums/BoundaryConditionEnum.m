classdef BoundaryConditionEnum

    properties (Constant = true)
    
        % Circular, wrapping, conditions over boundaries.
        CIRCULAR     = 'circ';       
        
        % Constant values along edges.
        CONST        = 'const';      
        
        % Fade arrays out to ICs at edges.
        FADE         = 'fade';       
        
        % Flips vector boundaries along shifting direction.
        FLIP         = 'flip';       
        
        % Linear interpolated boundary condition type.
        LINEAR       = 'linear';     
        
        % Mirrors across boundaries.
        MIRROR       = 'mirror';     
        
        % Transparent boundary condition type.
        TRANSPARENT  = 'trans';      
        
        % Immutable wall boundary set by initial conditions.
        WALL         = 'wall';       
        
        % Special BC mode for rotating axisymmetric objects.
        WORMHOLE     = 'wormhole';   
        
        % Zero fluxes but constant values.
        ZERO         = 'zero';
        
    end %PROPERTIES
end %CLASSDEF
