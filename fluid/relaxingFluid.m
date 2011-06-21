function relaxingFluid(run, mass, mom, ener, mag, grav, X)
%   This routine is responsible for the actual fluid fluxing. It utilizes the Relaxed, Total 
%   Variation Diminishing (TVD), Monotone Upwind Scheme for Fluid Dynamics (MUSCL) scheme for
%   non-linear fluxing terms. The relaxation technique ensures that the non-linear function
%   velocities are properly made upwind (See Jin & Xin '95, Ui-Li Pen '04).
%
%>< run      data manager object                                                      ImogenManager
%>< mass     mass density array (cell)                                                FluidArray
%>< mom      momentum density array (cell)                                            FluidArray(3)
%>< ener     energy density array (cell)                                              FluidArray
%>< mag      magnetic field (face)                                                    MagnetArray(3)
%>< grav     gravitational potential                                                  GravityArray
%>> X        vector index of current fluxing direction (1,2,or 3)                     int

    %--- Initialize ---%
    fluxFactor = run.time.dTime ./ run.DGRID{X};
    v          = [mass, mom(1), mom(2), mom(3), ener];

    YYY = 1;
    L = [X 2 3]; L(X)=1;
   
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%                   Half-Timestep predictor step (first-order upwind,not TVD)
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if run.useGPU == 1
        [freezea pressa] = freezeAndPtot(mass.array, ...
                                         ener.array, ...
                                         mom(L(1)).array, mom(L(2)).array, mom(L(3)).array, ...
                                         mag(L(1)).cellMag.array, mag(L(2)).cellMag.array, mag(L(3)).cellMag.array, ...
                                         run.GAMMA, 1);
        [v(1).store.array v(5).store.array ...
         v(L(1)+1).store.array v(L(2)+1).store.array v(L(3)+1).store.array] = cudaWstep(mass.array, ener.array, ...
                                                                         mom(L(1)).array, mom(L(2)).array, mom(L(3)).array, ...
                                                                         mag(L(1)).cellMag.array, mag(L(2)).cellMag.array, mag(L(3)).cellMag.array, ...
                                                                         pressa, freezea, fluxFactor, 1, run.pureHydro);

	cudaArrayAtomic(mass.store.array, run.fluid.MINMASS, ENUM.CUATOMIC_SETMIN);

else
    wFluidFlux(run, mass, mom, ener, mag, grav, run.fluid.freezeSpd(X), X);
    
    tempFreeze = run.fluid.freezeSpd(X).array;

    for i=1:5
        v(i).store.fluxR.array = 0.5*( v(i).array + v(i).wArray );
        v(i).store.fluxL.array = 0.5*( v(i).array - v(i).wArray );
        v(i).store.array = v(i).array - 0.5*fluxFactor .* tempFreeze .* ...
                         ( v(i).store.fluxR.array - v(i).store.fluxR.shift(YYY,-1) ...
                         + v(i).store.fluxL.array - v(i).store.fluxL.shift(YYY,1) );
        end
        v(i).store.cleanup();

    mass.store.array = max(mass.store.array, run.fluid.MINMASS);
end
   
%imagesc(double(v(1).store.array));
%input('cony: ');

    run.fluid.freezeSpd(X).cleanup();

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%                   Full-Timestep corrector step (second-order relaxed TVD)
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%    momStore = [mom(L(1)).store, mom(L(2)).store, mom(L(3)).store];
%    wFluidFlux(run, mass.store, momStore, ener.store, mag, grav, run.fluid.freezeSpdTVD(1), 1);
    
%    tempFreeze = run.fluid.freezeSpdTVD(1).array;

if run.useGPU

%1 = 1.0*mass.array;
%2 = 1.0*ener.array;
%3 = 1.0*mom(L(1)).array;
%r4 = 1.0*mom(L(2)).array;
%r5 = 1.0*mom(L(3)).array;


   [freezea pressa] = freezeAndPtot(mass.store.array, ener.store.array, ...
                                         mom(L(1)).store.array, mom(L(2)).store.array, mom(L(3)).store.array, ...
                                         mag(L(1)).cellMag.array, mag(L(2)).cellMag.array, mag(L(3)).cellMag.array, ...
                                         run.GAMMA, 1);

    cudaCorrectorStep(mass.store.array, ener.store.array, ...
                      mom(L(1)).store.array, mom(L(2)).store.array, mom(L(3)).store.array, ...
                      mag(L(1)).cellMag.array, mag(L(2)).cellMag.array, mag(L(3)).cellMag.array, ...
                      pressa, ...
                      mass.array, ener.array, mom(L(1)).array, mom(L(2)).array, mom(L(3)).array, ...
                      freezea, fluxFactor, 0);

    mass.applyStatics();
    ener.applyStatics();
    mom(1).applyStatics();
    mom(2).applyStatics();
    mom(3).applyStatics();
    
    cudaArrayAtomic(mass.array, run.fluid.MINMASS, ENUM.CUATOMIC_SETMIN);

return;
end

    for i=1:5
        if run.useGPU
            [v(i).fluxL.array v(i).fluxR.array] = cudaMHDKernels(8, v(i).store.array, v(i).store.wArray);
        else
        v(i).fluxR.array =  0.5*(v(i).store.array + v(i).store.wArray);
        v(i).fluxL.array =  0.5*(v(i).store.array - v(i).store.wArray);
        end
        
	% NOTE fluxes are not divided by two here, but in the flux limiters
        run.fluid.limiter{X}(v(i).fluxR, ...
                             (v(i).fluxR.array - v(i).fluxR.shift(YYY,-1)), ...
                             (v(i).fluxR.shift(YYY,1) - v(i).fluxR.array) );
        run.fluid.limiter{X}(v(i).fluxL, ...
                             (v(i).fluxL.shift(YYY,-1) - v(i).fluxL.array), ...
                             (v(i).fluxL.array - v(i).fluxL.shift(YYY,1)) );

        if run.useGPU
            v(i).array = cudaMHDKernels(7, v(i).array, tempFreeze, v(i).fluxR.array, v(i).fluxR.shift(YYY,-1), v(i).fluxL.array, v(i).fluxL.shift(YYY,1), fluxFactor);
        else
            v(i).array = v(i).array - fluxFactor .* tempFreeze .* ...
                                    ( v(i).fluxR.array  - v(i).fluxR.shift(YYY,-1) ...
                                    + v(i).fluxL.array  - v(i).fluxL.shift(YYY,1) );
        end

        v(i).cleanup();
    end

    if run.useGPU == 1
        cudaArrayAtomic(mass.store.array, run.fluid.MINMASS, ENUM.CUATOMIC_SETMIN);
    else
        mass.store.array = max(mass.store.array, run.fluid.MINMASS);
    end


    run.fluid.freezeSpdTVD(X).cleanup();

end
