function omegaAnalysis(analysis, fsize)

if nargin == 1; fsize = 16; end

vinv = 1 / analysis.equil.vel(1,1);

hold off; 
plot(vinv*analysis.omegaIm_fromdrho2./analysis.kyValues,'r');
hold on;
plot(vinv*analysis.omegaIm_fromdrho2./analysis.kyValues,'r-');
plot(vinv*analysis.omegaIm_fromdvx2 ./analysis.kyValues,'g-');
plot(vinv*analysis.omegaIm_fromdvy2 ./analysis.kyValues,'b-');
plot(vinv*analysis.omegaIm_fromdrho1./analysis.kyValues,'r-.');
plot(vinv*analysis.omegaIm_fromdvx1 ./analysis.kyValues,'g-.');
plot(vinv*analysis.omegaIm_fromdvy1 ./analysis.kyValues,'b-.');
plot(vinv*analysis.omegaIm_fromdby2 ./analysis.kyValues,'c-');
plot(vinv*analysis.omegaIm_fromdby1 ./analysis.kyValues,'c-.');
plot(vinv*imag(analysis.linear.omega(1:analysis.nModes(1)))'./analysis.kyValues,'k');
plot(0*numel(analysis.kyValues),'k-');

xlabel('Mode number','fontsize',fsize);
ylabel('imag(omega) / Ky Vin','fontsize',fsize);

title('Fits of omega* using multiple perturbed quantities','fontsize',fsize);

grid on

end
