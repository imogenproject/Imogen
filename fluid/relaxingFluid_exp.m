function relaxingFluid_exp(run, mass, mom, ener, mag, grav, X)
%   This routine is responsible for the actual fluid fluxing. It utilizes the Relaxed, Total 
%   Variation Diminishing (TVD), Monotone Upwind Scheme for Fluid Dynamics (MUSCL) scheme for
%   non-linear fluxing terms. The relaxation technique ensures that the non-linear function
%   velocities are properly made upwind (See Jin & Xin '95, Ui-Li Pen '04).
%
%>< run		data manager object														ImogenManager
%>< mass	mass density array (cell)												FluidArray
%>< mom		momentum density array (cell)											FluidArray(3)
%>< ener	energy density array (cell)												FluidArray
%>< mag		magnetic field (face)													MagnetArray(3)
%>< grav	gravitational potential													GravityArray
%>> X		vector index of current fluxing direction (1,2,or 3)					int

                        
    %-----------------------------------------------------------------------------------------------
    % Initialize
    %-----------
    fluxFactor = run.time.dTime ./ run.DGRID{X};
    v = [mass, mom(1), mom(2), mom(3), ener];
	
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%                   Half-Timestep predictor step (first-order upwind,not TVD)
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
	wFluidFlux(run, mass, mom, ener, mag, grav, run.fluid.freezeSpd(X), X);
	
    tempFreeze = run.fluid.freezeSpd(X).array;
	for i=1:5
		v(i).store.fluxR.array = 0.5*( v(i).array + v(i).wArray );
		v(i).store.fluxL.array = 0.5*( v(i).array - v(i).wArray );
		v(i).store.array = v(i).array - 0.5*fluxFactor .* tempFreeze .* ...
                        ( v(i).store.fluxR.array - v(i).store.fluxR.shift(X,-1) ...
                        + v(i).store.fluxL.array - v(i).store.fluxL.shift(X,1) );
		v(i).store.cleanup();
	end
	mass.store.array = max(mass.store.array, run.fluid.MINMASS);
	run.fluid.freezeSpd(X).cleanup();
		     
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%                   Full-Timestep corrector step (second-order relaxed TVD)
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	momStore = [mom(1).store, mom(2).store, mom(3).store];
	wFluidFlux(run, mass.store, momStore, ener.store, mag, grav, run.fluid.freezeSpdTVD(X), X);
    
    tempFreeze = run.fluid.freezeSpdTVD(X).array;
	for i=1:5
		v(i).fluxR.array =  0.5*(v(i).store.array + v(i).store.wArray);
		v(i).fluxR.array = vanleerLimiter( v(i).fluxR.array, ...
                                             0.5*(v(i).fluxR.array - v(i).fluxR.shift(X,-1)), ...
											 0.5*(v(i).fluxR.shift(X,1) - v(i).fluxR.array) );
                                        
    	v(i).fluxL.array =  0.5*(v(i).store.array - v(i).store.wArray);
		v(i).fluxL.array = vanleerLimiter( v(i).fluxL.array, ...
                                             0.5*(v(i).fluxL.shift(X,-1) - v(i).fluxL.array), ...
											 0.5*(v(i).fluxL.array - v(i).fluxL.shift(X,1)) );
                                        
    	v(i).array = v(i).array  - fluxFactor .* tempFreeze .* ...
                                ( v(i).fluxR.array  - v(i).fluxR.shift(X,-1) ...
                                + v(i).fluxL.array  - v(i).fluxL.shift(X,1) );
		v(i).cleanup();
	end

    mass.array = max(mass.array, run.fluid.MINMASS);
	run.fluid.freezeSpdTVD(X).cleanup();
	
end
