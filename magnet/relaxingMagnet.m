function relaxingMagnet(run, mag, velGrid, X, I)
% Flux the input magnetic field according to the two-step MUSCL TVD scheme and return the result as 
% well as the calculated flux to use later for the constraint quantity according to constrained 
% transport.
%
%><	run				run variable manager object								ImogenManager
%>< mag             magnetic field array (3D)								MagnetArray(3)
%>< velGrid         face aligned velocity array								FluidArray
%>> X				vector index for fluxing direction						int
%>>	I				the component of the b-field to operate on				int


    %-----------------------------------------------------------------------------------------------
    % Initialization
    %---------------
    fluxFactor = 0.5*run.time.dTime ./ run.DGRID{X};

    if run.useGPU
        [mag(I).store(X).array velocityFlow] = cudaMagW(mag(I).array, velGrid.array, fluxFactor, X);
    else

        velocityFlow = ( (velGrid.array + velGrid.shift(X,1)) < 0.0 );
        
        %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %Half-Timestep predictor step (first-order upwind,not TVD)
        %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
        mag(I).store(X).fluxR.array = mag(I).array .* velGrid.array;
        mag(I).store(X).fluxR.array = mag(I).store(X).fluxR.array .* (1-velocityFlow) ...
							  + mag(I).store(X).fluxR.shift(X,1) .* velocityFlow;
    
        mag(I).store(X).array = mag(I).array ...
				- fluxFactor .* (mag(I).store(X).fluxR.array - mag(I).store(X).fluxR.shift(X,-1));

    end

    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++
    %Full-Timestep corrector step (second-order relaxed TVD)
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    fluxFactor = 2*fluxFactor; %Multiply to get full timestep
    mag(I).wMag(X).array = mag(I).store(X).array .* velGrid.array;
    
    mag(I).flux(X).array = mag(I).wMag(X).array .* (1-velocityFlow) ...
                           + mag(I).wMag(X).shift(X,1) .* velocityFlow;
    dFluxR  = ( mag(I).wMag(X).shift(X,1) - mag(I).flux(X).array ) .* (1-velocityFlow) ...
            + ( mag(I).flux(X).array - mag(I).wMag(X).shift(X,2) ) .* velocityFlow;
    dFluxL  = ( mag(I).flux(X).array - mag(I).wMag(X).shift(X,-1) ) .* (1-velocityFlow) ...
            + ( mag(I).wMag(X).array - mag(I).flux(X).array ) .* velocityFlow;
    run.magnet.limiter{X}(mag(I).flux(X), dFluxL, dFluxR); % This is doubled, appropriate halving done by limiter functions

    mag(I).array = mag(I).array - fluxFactor .* ( mag(I).flux(X).array - mag(I).flux(X).shift(X,-1) );
    
    %-----------------------------------------------------------------------------------
    % Reuse advection flux for constraint step for CT
    %------------------------------------------------
    fluxFactor = run.time.dTime ./ run.DGRID{I};

    mag(I).flux(X).array = mag(I).flux(X).array - mag(I).flux(X).shift(I,1);
    mag(I).flux(X).array =  mag(I).flux(X).shift(X,-1);
    mag(X).array = mag(X).array - fluxFactor .* mag(I).flux(X).array;

end
