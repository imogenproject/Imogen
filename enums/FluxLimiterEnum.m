classdef FluxLimiterEnum
% Enumerated types for the available flux limiters that can be applied to the various flux routines.

    properties (Constant = true)
    
        % Van Leer flux limiter type.
        VAN_LEER     = 'vanleer';
        
        % Superbee flux limiter type.
        SUPERBEE     = 'superbee';
        
        % MinMod flux limiter type.
        MINMOD       = 'minmond';    
    
    end %PROPERTIES


end %CLASSDEF
