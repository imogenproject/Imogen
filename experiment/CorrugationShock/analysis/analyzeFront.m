function [growthrates growresidual phaserates phaseresidual] = analyzeFront(frontFFT, timeVals, linearFrames);  
% >> frontFFT      Fourier transform of the front's position in space
% >> timeVals      The simulation time at each step saved
% >> linearFrames  The set of saved timesteps that the simulation was found to be linear in
%
% << growthrates   The calculated growth rates for the front perturbation
% << correlation   The residual norm for the growth rate fit (lower = better)

% Iterate over modes
for u = 1:size(frontFFT, 1); for v = 1:size(frontFFT,2)
	% Fit a linear polynomial to the shock perturbation's amplitude for that mode.
        [f s]= polyfit(squeeze(timeVals(linearFrames)), squeeze(log(abs(frontFFT(u,v,linearFrames))))', 1);
        growthrates(u,v) = f(1);
        growresidual(u,v) = s.normr;

        [f s]= polyfit(squeeze(timeVals(linearFrames)), unwrap(squeeze(angle(frontFFT(u,v,linearFrames)))'), 1);
        phaserates(u,v) = f(1);
        phaseresidual(u,v) = s.normr;
      
end; end

end
