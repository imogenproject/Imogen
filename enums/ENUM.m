classdef ENUM
    
    properties (Constant = true)
        MASS      = 'mass';
        MOM       = 'mom';
        ENER      = 'ener';
        MAG       = 'mag';
        GRAV      = 'grav';
        
        SCALAR    = 0;
        VECTOR    = [1 2 3];
        
        POINT_FADE = 'point';
        
        %--- Pressure Types ---%
        PRESSURE_TOTAL_AND_SOUND = 'totsnd';
        PRESSURE_SOUND_SPEED     = 'sound';
        PRESSURE_GAS             = 'gas'
        PRESSURE_TOTAL           = 'total'
        PRESSURE_MAGNETIC        = 'magnetic'
        
        %--- Boundary Condition Modes  ---%
        BCMODE_CIRCULAR     = 'circ';       % Circular, wrapping, conditions over boundaries.
        BCMODE_CONST        = 'const';      % Constant values along edges.
        BCMODE_FADE         = 'fade';       % Fade arrays out to ICs at edges.
        BCMODE_FLIP         = 'flip';       % Flips vector boundaries along shifting direction.
        BCMODE_LINEAR       = 'linear';     % Linear interpolated boundary condition type.
        BCMODE_MIRROR       = 'mirror';     % Mirrors across boundaries.
        BCMODE_TRANSPARENT  = 'trans';      % Transparent boundary condition type.
        BCMODE_WALL         = 'wall';       % Immutable wall boundary set by initial conditions.
        BCMODE_WORMHOLE     = 'wormhole';   % Special BC mode for rotating axisymmetric objects.
        BCMODE_ZERO         = 'zero';       % Zero fluxes but constant values.

        %--- Time update Modes ---%
        TIMEUPDATE_PER_STEP       = 0;      % Update timestep every flux step.
        TIMEUPDATE_PER_ITERATION  = 1;      % Update timestemp every iteration.

        %--- Gravitational Solvers ---%
        GRAV_SOLVER_EMPTY       = 'empty';         % Empty solver.
        GRAV_SOLVER_BICONJ      = 'biconj';        % Linear solver using BICGSTAB.
        GRAV_SOLVER_GPU         = 'biconj_gpu';    % Linear solver using BICGSTAB running on GPU
        GRAV_SOLVER_MULTIGRID   = 'multigrid';     % Hierarchial discretization solver.
        
        GRAV_BCSOURCE_FULL      = 'full';          % Computes boundary conditions at every boundary cell
        GRAV_BCSOURCE_INTERP    = 'interpolated';  % Computes boundary conditions at every 4th boundary
                                                   % cell and interpolates

        %--- Radiation Model Types ---%
        RADIATION_NONE              = 'empty';
        RADIATION_OPTICALLY_THIN    = 'optically_thin';
        
        %--- Artificial Viscosity Types ---%
        ARTIFICIAL_VISCOSITY_NONE               = 'empty';
        ARTIFICIAL_VISCOSITY_NEUMANN_RICHTMYER  = 'neumann_richtmyer';
        ARTIFICIAL_VISCOSITY_CARAMANA_SHASHKOV_WHALEN = 'caramana_shashkov_whalen';

	CUATOMIC_SETMIN  = 1;
	CUATOMIC_SETMAX  = 2;
	CUATOMIC_FIXNAN  = 3;
        
    end

end
