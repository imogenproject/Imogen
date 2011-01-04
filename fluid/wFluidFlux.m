function wFluidFlux(run, mass, mom, ener, mag, grav, freezeSpd, X)
% This function populates the auxiliary, wArray, arrays for the fluid variables as needed for the
% fluid fluxing routine.
%
%>< run             data manager object                                            ImogenManager
%>< mass            mass density array (cell)                                      FluidArray
%>< mom             momentum density array (vector,cell)                           FluidArray(3)
%>< ener            energy density array (cell)                                    FluidArray
%>< mag             magnetic field                                                 MagnetArray(3)
%>< grav            gravitational potential                                        GravityArray
%>< freezeSpd       freeze speed object                                            InitializedArray
%>> X               vector index of current fluxing direction (1,2,or 3)           int

    dirVec    = zeros(3,1); dirVec(X)  = 1.0;
    repsMat = ones(3,1);  repsMat(X) = run.gridSize(X);

    %--- Determine total pressure and relaxed freezing speed ---%
    [press, soundSpd] = pressure('totsnd',run, mass, mom, ener, mag);
    press             = press + run.fluid.viscosity.solve(run, mass, mom, ener, mag, soundSpd, X);  % (run, mass, mom, soundSpd, X);
    velocity          = mom(X).array ./ mass.array;


    if run.useGPU
        spd               = GPUdouble(max( double(abs(velocity) + soundSpd) ,[],X) );
    else
        spd               = max( (abs(velocity) + soundSpd) ,[],X);
    end

    if iscodistributed(spd), spd = gather(spd); end

    freezeSpd.array = repmat( spd, repsMat' );

    %-----------------------------------------------------------------------------------------------
    % Determine the auxiliary relaxing function arrays
    %-------------------------------------------------
    
    %--- MASS DENSITY ---%
    mass.wArray    = mom(X).array ./ freezeSpd.array;
    


    %--- ENERGY DENSITY ---%
    ener.wArray    = velocity .* (ener.array + press) - mag(X).cellMag.array .* ...
                        ( mag(1).cellMag.array .* mom(1).array ...
                        + mag(2).cellMag.array .* mom(2).array ...
                        + mag(3).cellMag.array .* mom(3).array) ./ mass.array;
    ener.wArray    = ener.wArray ./ freezeSpd.array;
    
    %--- MOMENTUM DENSITY ---%
    for i=1:3
        mom(i).wArray    = velocity .* mom(i).array + press*dirVec(i)...
                             - mag(X).cellMag.array .* mag(i).cellMag.array;
        mom(i).wArray    = mom(i).wArray ./ freezeSpd.array;

q=isnan(double(mom(i).wArray));
if max(q(:))>0;
mass.array
velocity
%mom(i).array
%press

error('mom nan'); end
    end

end
