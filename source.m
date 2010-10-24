function source(run, mass, mom, ener, mag, grav)
% This function sources the non-conservative terms in the MHD equations like gravitational potential
% and radiation terms. Effectively it provides the means to add terms that cannot be brought within 
% the del operator, which is the foundation of the spatial fluxing routines.
%
%>< run			data manager object.                                            ImogenManager
%>< mass		mass density                                                    FluidArray  
%>< mom			momentum density                                                FluidArray(3)
%>< ener        energy density                                                  FluidArray
%>< mag         magnetic field density                                          FluidArray(3)
%>< grav		gravitational potential                                         GravityArray

    %--- Gravitational Potential Sourcing ---%
    %       If the gravitational portion of the code is active, the gravitational potential terms
    %       in both the momentum and energy equations must be appended as source terms.
    if run.gravity.ACTIVE
        enerSource = zeros(run.gridSize);
        for i=1:3
            momSource       = run.time.dTime*mass.thresholdArray ...
                                                    .* grav.calculate5PtDerivative(i,run.DGRID{i});
            enerSource      = enerSource + momSource .* mom(i).array ./ mass.array;
            mom(i).array    = mom(i).array - momSource;
        end
        ener.array          = ener.array - enerSource;
    end
    
    %--- Radiation Sourcing ---%
    %       If radiation is active, the radiation terms are subtracted, as a sink, from the energy
    %       equation.
    ener.array              = ener.array - run.fluid.radiation.solve(run, mass, mom, ener, mag);
    
end
