function plotFrontAnalysis_2d(front, timeVals, linearFrames, fixfsize)
% Generates a nice plot using the front tracking information outputted by the growth analyzer.


amp   = squeeze(log(abs(front.FFT)));
phase = squeeze(angle(front.FFT));

wavenum = (1:size(amp,1))-1;
wavevec = wavenum * 2*pi/(.01*size(amp,1));
fsize=16; if nargin == 4; fsize = fixfsize; end
figure();

reliableWaves = 1:round(numel(wavenum)/10);
kMax = round(numel(wavenum)/10);

subpltSize = .35;
subplotOffset = .07;

subplot(2,2,1);
%size(front.growthRate(reliableWaves))
%size(wavenum(reliableWaves))
	plot(wavenum(reliableWaves), front.growthRate(reliableWaves)./wavenum(reliableWaves), 'b');

	sFactor = max(front.growthRate(reliableWaves)) / max(front.residualNorm(reliableWaves));
	hold on;
	plot(wavenum(reliableWaves), front.residualNorm(reliableWaves)*sFactor*.75, 'r');
	hold off;
	grid on;

	xlabel('Wavevector','fontsize',fsize);
	ylabel('B: Growth rate linear fit; R: residual norm','fontsize',fsize);
	title('Growth rate and error norms','fontsize',fsize);
	set(gca,'position',[subplotOffset .5+subplotOffset subpltSize subpltSize]);

subplot(2,2,2);
	a = round(kMax/3);
	b = round(2*kMax/3);

	plot(timeVals, amp(2:a,1:end),'r');
        hold on;
        plot(timeVals, amp((a+1):b,1:end),'g');
        plot(timeVals, amp((b+1):kMax,1:end),'b');
	hold off;

	xlabel('Time, simulation units');
	ylabel('log(abs(fourier amplitude))','fontsize',fsize);

	title('Fourier spectrum of shock surface''s position over time','fontsize',fsize);
	set(gca,'position', [.5+subplotOffset .5+subplotOffset subpltSize subpltSize]);

subplot(2,2,3);
	imagesc(phase(reliableWaves,:));
    labelstring = sprintf('Frame #; Times 0 to %.4g', timeVals(end));
    xlabel(labelstring,'fontsize',fsize);
	ylabel('Mode # + 1','fontsize',fsize);
	title('Phase of shock front''s modes','fontsize',fsize);
	set(gca,'position', [subplotOffset subplotOffset subpltSize subpltSize]);

subplot(2,2,4);
	plot(timeVals, unwrap(phase(round(kMax/3):round(kMax*2/3),1:end)' ));
	
	grid on;
	xlabel('Time, simulation units','fontsize',fsize);
	ylabel('Mode phase, radians','fontsize',fsize);
	title(['Phase over time of modes 1-' sprintf('%i',kMax)],'fontsize',fsize);
	set(gca,'position', [.5+subplotOffset subplotOffset subpltSize subpltSize]);


end 
