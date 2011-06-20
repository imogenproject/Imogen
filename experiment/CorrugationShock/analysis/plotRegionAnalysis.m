function plotRegionAnalysis(front, damping, timeVals, lastframe, fixfsize)
% Generates a nice plot using the front tracking information outputted by the growth analyzer.

figure();
fsize=16; if nargin == 5; fsize = fixfsize; end

subpltSize = .35;
subplotOffset = .07;

subplot(2,2,1);
	plot(timeVals, log(front.rms));

	xlabel('Time, simulation units','fontsize',fsize);
	ylabel('RMS front distortion','fontsize',fsize);
	title('RMS shock front perturbation amplitude','fontsize',fsize);
	set(gca,'position', [subplotOffset .5+subplotOffset subpltSize subpltSize]);

subplot(2,2,2);
    waveVec = damping.KY;

    plot(waveVec, damping.dampPre./waveVec,'g');
    hold on;
sfact = max(max(abs(damping.corrPre)), max(abs(damping.corrPost))) / ...
        max(max(abs(damping.dampPre./waveVec)), max(abs(damping.dampPost./waveVec)));
    
    plot(waveVec, damping.corrPre/sfact,'g.');
    plot(waveVec, damping.dampPost./waveVec,'b')
    plot(waveVec, damping.corrPost/sfact,'b.');
    
	plot(waveVec, zeros(size(waveVec)),'k-');
    
	xlabel('Wavevector','fontsize',fsize);
	ylabel('B: im(kxpost*), G: im(kxpre*)','fontsize',fsize);

	title('Ky-normalized imaginary part of postshock waves','fontsize',fsize);
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
