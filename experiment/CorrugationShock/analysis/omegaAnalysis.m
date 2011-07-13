function omegaAnalysis(analysis, fsize)

if nargin == 1; fsize = 16; end

vinv = 1 / analysis.equil.vel(1,1);

hold off; 

subplot(1,2,1);

plot(vinv*imag(analysis.omega_fromdrho2)./analysis.kyWavenums,'r','DisplayName','drho post');
hold on;
plot(vinv*imag(analysis.omega_fromdrho2)./analysis.kyWavenums,'r-','DisplayName','drho post');
plot(vinv*imag(analysis.omega_fromdvx2) ./analysis.kyWavenums,'g-','DisplayName','dvx post');
plot(vinv*imag(analysis.omega_fromdvy2) ./analysis.kyWavenums,'b-','DisplayName','dvy post');
plot(vinv*imag(analysis.omega_fromdrho1)./analysis.kyWavenums,'r-.','DisplayName','drho pre');
plot(vinv*imag(analysis.omega_fromdvx1) ./analysis.kyWavenums,'g-.','DisplayName','dvx pre');
plot(vinv*imag(analysis.omega_fromdvy1) ./analysis.kyWavenums,'b-.','DisplayName','dvy pre');
plot(vinv*imag(analysis.omega_fromdby2) ./analysis.kyWavenums,'c-','DisplayName','dby post');
plot(vinv*imag(analysis.omega_fromdby1) ./analysis.kyWavenums,'c-.','DisplayName','dby pre');
plot(vinv*imag(analysis.linear.omega(1:analysis.nModes(1)))./analysis.kyWavenums','k','DisplayName','eta');
plot(0*numel(analysis.kyWavenums),'k-');

xlabel('Mode number','fontsize',fsize);
ylabel('imag(omega) / Ky Vin','fontsize',fsize);

title('Fits of im(omega*) using multiple perturbed quantities','fontsize',fsize);

grid on

subplot(1,2,2);

plot(vinv*real(analysis.omega_fromdrho2)./analysis.kyWavenums,'r','DisplayName','drho post');
hold on;
plot(vinv*real(analysis.omega_fromdrho2)./analysis.kyWavenums,'r-','DisplayName','drho post');
plot(vinv*real(analysis.omega_fromdvx2) ./analysis.kyWavenums,'g-','DisplayName','dvx post');
plot(vinv*real(analysis.omega_fromdvy2) ./analysis.kyWavenums,'b-','DisplayName','dvy post');
plot(vinv*real(analysis.omega_fromdrho1)./analysis.kyWavenums,'r-.','DisplayName','drho pre');
plot(vinv*real(analysis.omega_fromdvx1) ./analysis.kyWavenums,'g-.','DisplayName','dvx pre');
plot(vinv*real(analysis.omega_fromdvy1) ./analysis.kyWavenums,'b-.','DisplayName','dvy pre');
plot(vinv*real(analysis.omega_fromdby2) ./analysis.kyWavenums,'c-','DisplayName','dby post');
plot(vinv*real(analysis.omega_fromdby1) ./analysis.kyWavenums,'c-.','DisplayName','dby pre');
plot(vinv*real(analysis.linear.omega(1:analysis.nModes(1)))./analysis.kyWavenums','k');
plot(0*numel(analysis.kyWavenums),'k-');

grid on;

xlabel('Mode number','fontsize',fsize);
ylabel('real(omega) / Ky Vin','fontsize',fsize);

title('Fits of re(omega*) using multiple perturbed quantities','fontsize',fsize);


end
