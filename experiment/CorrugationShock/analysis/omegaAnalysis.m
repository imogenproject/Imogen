function omegaAnalysis(analysis, fsize)

if nargin == 1; fsize = 16; end

vinv = 1 / norm(analysis.equil.vel(:,1));

hold off; 
subplot(1,2,1);

plot(vinv*imag(analysis.omega.fromdrho2)./analysis.kyValues,'r','DisplayName','drho post');
hold on;
plot(vinv*imag(analysis.omega.fromdrho2)./analysis.kyValues,'r-','DisplayName','drho post');
plot(vinv*imag(analysis.omega.fromdvx2) ./analysis.kyValues,'g-','DisplayName','dvx post');
plot(vinv*imag(analysis.omega.fromdvy2) ./analysis.kyValues,'b-','DisplayName','dvy post');
plot(vinv*imag(analysis.omega.fromdbx2) ./analysis.kyValues,'m-','DisplayName','dbx post');
plot(vinv*imag(analysis.omega.fromdby2) ./analysis.kyValues,'c-','DisplayName','dby post');

plot(vinv*imag(analysis.omega.fromdrho1)./analysis.kyValues,'r-.','DisplayName','drho pre');
plot(vinv*imag(analysis.omega.fromdvx1) ./analysis.kyValues,'g-.','DisplayName','dvx pre');
plot(vinv*imag(analysis.omega.fromdvy1) ./analysis.kyValues,'b-.','DisplayName','dvy pre');
plot(vinv*imag(analysis.omega.fromdbx1) ./analysis.kyValues,'m-.','DisplayName','dbx pre')
plot(vinv*imag(analysis.omega.fromdby1) ./analysis.kyValues,'c-.','DisplayName','dby pre');

plot(vinv*imag(analysis.omega.front(1:analysis.nModes(1)))./analysis.kyValues','k','DisplayName','eta');

plot(1:analysis.nModes(1), zeros(size(analysis.kyValues)),'k-','DisplayName','');

xlabel('Mode number','fontsize',fsize);
ylabel('imag(omega) / Ky Vin','fontsize',fsize);

title('Fits of im(omega*) using multiple perturbed quantities','fontsize',fsize);

grid on

subplot(1,2,2);
plot(vinv*real(analysis.omega.fromdrho2)./analysis.kyValues,'r','DisplayName','drho post');
hold on;
plot(vinv*real(analysis.omega.fromdrho2)./analysis.kyValues,'r-','DisplayName','drho post');
plot(vinv*real(analysis.omega.fromdvx2) ./analysis.kyValues,'g-','DisplayName','dvx post');
plot(vinv*real(analysis.omega.fromdvy2) ./analysis.kyValues,'b-','DisplayName','dvy post');
plot(vinv*real(analysis.omega.fromdbx2) ./analysis.kyValues,'m-','DisplayName','dbx post');
plot(vinv*real(analysis.omega.fromdby2) ./analysis.kyValues,'c-','DisplayName','dby post');

plot(vinv*real(analysis.omega.fromdrho1)./analysis.kyValues,'r-.','DisplayName','drho pre');
plot(vinv*real(analysis.omega.fromdvx1) ./analysis.kyValues,'g-.','DisplayName','dvx pre');
plot(vinv*real(analysis.omega.fromdvy1) ./analysis.kyValues,'b-.','DisplayName','dvy pre');
plot(vinv*real(analysis.omega.fromdbx1) ./analysis.kyValues,'m-.','DisplayName','dbx pre')
plot(vinv*real(analysis.omega.fromdby1) ./analysis.kyValues,'c-.','DisplayName','dby pre');

plot(vinv*real(analysis.omega.front(1:analysis.nModes(1)))./analysis.kyValues','k','DisplayName','eta');

plot(1:analysis.nModes(1), zeros(size(analysis.kyValues)),'k-','DisplayName','');

grid on;

xlabel('Mode number','fontsize',fsize);
ylabel('real(omega) / Ky Vin','fontsize',fsize);

title('Fits of re(omega*) using multiple perturbed quantities','fontsize',fsize);


end
