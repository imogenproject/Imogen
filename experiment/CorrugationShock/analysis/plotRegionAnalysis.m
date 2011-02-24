function plotRegionAnalysis(front, searchRegion, timeVals, lastframe, linearFrames)
% Generates a nice plot using the front tracking information outputted by the growth analyzer.

figure();

subplot(2,2,1);
	plot(timeVals, log(front.rms));

	xlabel('Time, simulation units');
	ylabel('RMS front distortion');
	title('RMS shock front perturbation amplitude');
	set(gca,'position', [.05 .55 .4 .4]);

subplot(2,2,2);
	plot(timeVals, searchRegion.mdr,'b');

	dev = searchRegion.corr; dev(isnan(dev)) = 0;
	dev = sum(dev');
	sFactor = max(abs(dev)) / max(searchRegion.mdr);
	hold on;
	plot(timeVals, dev / sFactor,'r');
	plot(timeVals, zeros(size(timeVals)),'k');

	xlabel('Time, simulation units');
	ylabel('B: Exponent; Red: residual');

	title('Fitted decay rate of selected region perturbations');
	set(gca,'position', [.55 .55 .4 .4]);

subplot(2,2,3);
	x = size(lastframe.mass,1);
	y = size(lastframe.mass,2);

	surf(lastframe.mass(round(x/2+6):min(round(x/2+x/6+6),round(x/2)+100), 1:min(256, y), 1)-lastframe.mass(end,1,1),'linestyle','none');
	title('Postshock density perturbation');
	set(gca,'position', [.05 .05 .4 .4]);

subplot(2,2,4);
        x = size(lastframe.mass,1);
        y = size(lastframe.mass,2);

        surf(lastframe.mass(max(round(x/2-x/6-6),round(x/2)-100):round(x/2-6), 1:min(256, y), 1)-lastframe.mass(1,1,1),'linestyle','none');
        title('Preshock density perturbation');
	set(gca,'position', [.55 .05 .4 .4]);


end 
