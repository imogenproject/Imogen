classdef PressureTypesEnum
% Pressure type enumerations passed as an argument to the  pressure routine to specify what type 
% of pressure to return.


    properties (Constant = true)
    
        % Enumerated type specifying both total pressure and sound speed.
        TOTAL_AND_SOUND = 'totsnd';
        
        % Enumerated type specifying the sound pressure.
        SOUND_SPEED     = 'sound';
        
        % Enumerated type specifying the gas pressure.
        GAS             = 'gas';
        
        % Enumerated type specifying the total pressure (gas + magnetic).
        TOTAL           = 'total';
        
        % Enumerated type specifying the magnetic pressure.
        MAGNETIC        = 'magnetic';
    
    
    end %PROPERTIES


end %CLASSDEF
