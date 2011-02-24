function plotFrontAnalysis_2d(front, timeVals, linearFrames)
% Generates a nice plot using the front tracking information outputted by the growth analyzer.


amp   = squeeze(log(abs(front.FFT)));
phase = squeeze(angle(front.FFT));

wavenum = 1:size(amp,1);
wavevec = wavenum * 2*pi/(.01*size(amp,1));

figure();

reliableWaves = 1:round(numel(wavenum)/15);
kMax = round(numel(wavenum)/15);

subplot(2,2,1);
	plot(wavenum(reliableWaves)-1, front.growthRate(reliableWaves), 'b');

	sFactor = max(front.growthRate(reliableWaves)) / max(front.residualNorm(reliableWaves));
	hold on;
	plot(wavenum(reliableWaves)-1, front.residualNorm(reliableWaves)*sFactor*.75, 'r');
	hold off;
	grid on;

	xlabel('Wavevector');
	ylabel('B: Growth rate linear fit; R: residual norm');
	title('Growth rate and error norms');
	set(gca,'position', [.05 .55 .4 .4]);

subplot(2,2,2);
	a = round(kMax/3);
	b = round(2*kMax/3);

	plot(timeVals, amp(2:a,2:end),'r');
        hold on;
        plot(timeVals, amp((a+1):b,2:end),'g');
        plot(timeVals, amp((b+1):kMax,2:end),'b');
	hold off;

	xlabel('Time, simulation units');
	ylabel('log(abs(fourier amplitude))');

	title('Fourier spectrum of shock surface''s position over time');
	set(gca,'position', [.55 .55 .4 .4]);

subplot(2,2,3);
	imagesc(phase(reliableWaves,:));
	xlabel(['Frame no; timespan 0 to ' timeVals(end)]);
	ylabel('Mode # + 1');
	title('Phase of shock front''s modes');
	set(gca,'position', [.05 .05 .4 .4]);

subplot(2,2,4);
	plot(timeVals, phase(round(kMax/3):round(kMax*2/3),2:end)');
	
	grid on;
	xlabel('Time, simulation units');
	ylabel('Mode phase, radians');
	title(['Phase over time of modes 1-' kMax]);
	set(gca,'position', [.55 .05 .4 .4]);


end 
