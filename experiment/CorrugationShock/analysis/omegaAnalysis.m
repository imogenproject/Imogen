function omegaAnalysis(analysis, fsize)

if nargin == 1; fsize = 16; end

vinv = 1 / analysis.equil.vel(1,1);

hold off; 

subplot(1,2,1);

plot(vinv*imag(analysis.omega_fromdrho2)./analysis.kyValues,'r');
hold on;
plot(vinv*imag(analysis.omega_fromdrho2)./analysis.kyValues,'r-');
plot(vinv*imag(analysis.omega_fromdvx2) ./analysis.kyValues,'g-');
plot(vinv*imag(analysis.omega_fromdvy2) ./analysis.kyValues,'b-');
plot(vinv*imag(analysis.omega_fromdrho1)./analysis.kyValues,'r-.');
plot(vinv*imag(analysis.omega_fromdvx1) ./analysis.kyValues,'g-.');
plot(vinv*imag(analysis.omega_fromdvy1) ./analysis.kyValues,'b-.');
plot(vinv*imag(analysis.omega_fromdby2) ./analysis.kyValues,'c-');
plot(vinv*imag(analysis.omega_fromdby1) ./analysis.kyValues,'c-.');
plot(vinv*imag(analysis.linear.omega(1:analysis.nModes(1)))./analysis.kyValues','k');
plot(0*numel(analysis.kyValues),'k-');

xlabel('Mode number','fontsize',fsize);
ylabel('imag(omega) / Ky Vin','fontsize',fsize);

title('Fits of im(omega*) using multiple perturbed quantities','fontsize',fsize);

grid on

subplot(1,2,2);

plot(vinv*real(analysis.omega_fromdrho2)./analysis.kyValues,'r');
hold on;
plot(vinv*real(analysis.omega_fromdrho2)./analysis.kyValues,'r-');
plot(vinv*real(analysis.omega_fromdvx2) ./analysis.kyValues,'g-');
plot(vinv*real(analysis.omega_fromdvy2) ./analysis.kyValues,'b-');
plot(vinv*real(analysis.omega_fromdrho1)./analysis.kyValues,'r-.');
plot(vinv*real(analysis.omega_fromdvx1) ./analysis.kyValues,'g-.');
plot(vinv*real(analysis.omega_fromdvy1) ./analysis.kyValues,'b-.');
plot(vinv*real(analysis.omega_fromdby2) ./analysis.kyValues,'c-');
plot(vinv*real(analysis.omega_fromdby1) ./analysis.kyValues,'c-.');
plot(vinv*real(analysis.linear.omega(1:analysis.nModes(1)))./analysis.kyValues','k');
plot(0*numel(analysis.kyValues),'k-');

grid on;

xlabel('Mode number','fontsize',fsize);
ylabel('real(omega) / Ky Vin','fontsize',fsize);

title('Fits of re(omega*) using multiple perturbed quantities','fontsize',fsize);


end
