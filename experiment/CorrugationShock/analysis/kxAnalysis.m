function omegaAnalysis(analysis, fsize)

if nargin == 1; fsize = 16; end

vinv = 1 / analysis.equil.vel(1,1);

hold off; 

subplot(1,2,1);

plot(vinv*imag(analysis.post.drhoKx)./analysis.kyValues,'r');
hold on;
plot(vinv*imag(analysis.pre.drhoKx)./analysis.kyValues,'r-');
plot(vinv*imag(analysis.post.dvxKx) ./analysis.kyValues,'g-');
plot(vinv*imag(analysis.post.dvyKx) ./analysis.kyValues,'b-');
plot(vinv*imag(analysis.pre.drhoKx)./analysis.kyValues,'r-.');
plot(vinv*imag(analysis.pre.dvxKx) ./analysis.kyValues,'g-.');
plot(vinv*imag(analysis.pre.dvyKx) ./analysis.kyValues,'b-.');
plot(vinv*imag(analysis.pre.dbyKx) ./analysis.kyValues,'c-');
plot(vinv*imag(analysis.pre.dbyKx) ./analysis.kyValues,'c-.');
plot(vinv*real(analysis.linear.omega(1:analysis.nModes(1)))./analysis.kyValues','k');
plot(0*numel(analysis.kyValues),'k-');

xlabel('Mode number','fontsize',fsize);
ylabel('imag(kx) / Ky','fontsize',fsize);

title('Fits of im(kx/ky) using multiple perturbed quantities','fontsize',fsize);

grid on

subplot(1,2,2);

plot(vinv*real(analysis.post.drhoKx)./analysis.kyValues,'r');
hold on;
plot(vinv*real(analysis.pre.drhoKx)./analysis.kyValues,'r-');
plot(vinv*real(analysis.post.dvxKx) ./analysis.kyValues,'g-');
plot(vinv*real(analysis.post.dvyKx) ./analysis.kyValues,'b-');
plot(vinv*real(analysis.pre.drhoKx)./analysis.kyValues,'r-.');
plot(vinv*real(analysis.pre.dvxKx) ./analysis.kyValues,'g-.');
plot(vinv*real(analysis.pre.dvyKx) ./analysis.kyValues,'b-.');
plot(vinv*real(analysis.pre.dbyKx) ./analysis.kyValues,'c-');
plot(vinv*real(analysis.pre.dbyKx) ./analysis.kyValues,'c-.');
plot(vinv*real(analysis.linear.omega(1:analysis.nModes(1)))./analysis.kyValues','k');
plot(0*numel(analysis.kyValues),'k-');

grid on;

xlabel('Mode number','fontsize',fsize);
ylabel('real(kx) / Ky','fontsize',fsize);

title('Fits of re(kx/ky) using multiple perturbed quantities','fontsize',fsize);


end
