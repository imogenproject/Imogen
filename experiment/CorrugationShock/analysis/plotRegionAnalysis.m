function plotRegionAnalysis(front, searchRegion, timeVals, lastframe, linearFrames, fixfsize)
% Generates a nice plot using the front tracking information outputted by the growth analyzer.

figure();
fsize=16; if nargin == 6; fsize = fixfsize; end

subpltSize = .35;
subplotOffset = .07;

subplot(2,2,1);
	plot(timeVals, log(front.rms));

	xlabel('Time, simulation units','fontsize',fsize);
	ylabel('RMS front distortion','fontsize',fsize);
	title('RMS shock front perturbation amplitude','fontsize',fsize);
	set(gca,'position', [subplotOffset .5+subplotOffset subpltSize subpltSize]);

subplot(2,2,2);
	plot(timeVals, searchRegion.mdrPre,'r');
	dev = searchRegion.corrPre; dev(isnan(dev)) = 0;
	dev = sum(dev');
	sFactor = max(abs(dev)) / max(searchRegion.mdrPre);
	hold on;
	plot(timeVals, dev / sFactor,'r-.');
   
    plot(timeVals, searchRegion.mdrPost,'g');
	dev = searchRegion.corrPost; dev(isnan(dev)) = 0;
	dev = sum(dev');
	sFactor = max(abs(dev)) / max(searchRegion.mdrPost);
	hold on;
	plot(timeVals, dev / sFactor,'g-.');
    
	plot(timeVals, zeros(size(timeVals)),'k');
    
	xlabel('Time, simulation units','fontsize',fsize);
	ylabel('B: Exponent; Red: residual','fontsize',fsize);

	title('Fitted decay rate of pre (r) and post(g) waves','fontsize',fsize);
	set(gca,'position', [.5+subplotOffset .5+subplotOffset subpltSize subpltSize]);

subplot(2,2,3);
	x = size(lastframe.mass,1);
	y = size(lastframe.mass,2);

	surf(lastframe.mass(round(x/2+6):min(round(x/2+x/6+6),round(x/2)+100), 1:min(256, y), 1)-lastframe.mass(end,1,1),'linestyle','none');
	title('Postshock density perturbation','fontsize',fsize);
	set(gca,'position', [subplotOffset subplotOffset subpltSize subpltSize]);

subplot(2,2,4);
        x = size(lastframe.mass,1);
        y = size(lastframe.mass,2);

        surf(lastframe.mass(max(round(x/2-x/6-6),round(x/2)-100):round(x/2-6), 1:min(256, y), 1)-lastframe.mass(1,1,1),'linestyle','none');
        title('Preshock density perturbation','fontsize',fsize);
	set(gca,'position', [.5+subplotOffset subplotOffset subpltSize subpltSize]);


end 
