function omegaAnalysis(analysis, fsize)

if nargin == 1; fsize = 16; end

hold off; 

subplot(1,2,1);

plot(imag(analysis.post.drhoKx)./analysis.kyWavenums,'r','DisplayName','drho post');
hold on;
plot(imag(analysis.post.dvxKx) ./analysis.kyWavenums,'g-','DisplayName','dvx post');
plot(imag(analysis.post.dvyKx) ./analysis.kyWavenums,'b-','DisplayName','dvy post');
plot(imag(analysis.post.dbyKx) ./analysis.kyWavenums,'c-','DisplayName','dby post');

plot(imag(analysis.pre.drhoKx)./analysis.kyWavenums,'r-.','DisplayName','drho pre');
plot(imag(analysis.pre.dvxKx) ./analysis.kyWavenums,'g-.','DisplayName','dvx pre');
plot(imag(analysis.pre.dvyKx) ./analysis.kyWavenums,'b-.','DisplayName','dvy pre');
plot(imag(analysis.pre.dbyKx) ./analysis.kyWavenums,'c-.','DisplayName','dby pre');
plot(real(analysis.linear.omega(1:analysis.nModes(1)))./analysis.kyWavenums','k');
plot(0*numel(analysis.kyWavenums),'k-');

xlabel('Mode number','fontsize',fsize);
ylabel('imag(kx) / Ky','fontsize',fsize);

title('Fits of im(kx/ky) using multiple perturbed quantities','fontsize',fsize);

grid on

subplot(1,2,2);

plot(real(analysis.post.drhoKx)./analysis.kyWavenums,'r','DisplayName','drho post');
hold on;
plot(real(analysis.post.dvxKx) ./analysis.kyWavenums,'g-','DisplayName','dvx post');
plot(real(analysis.post.dvyKx) ./analysis.kyWavenums,'b-','DisplayName','dvy post');
plot(real(analysis.post.dbyKx) ./analysis.kyWavenums,'c-','DisplayName','dby post');

plot(real(analysis.pre.drhoKx)./analysis.kyWavenums,'r-.','DisplayName','drho pre');
plot(real(analysis.pre.dvxKx) ./analysis.kyWavenums,'g-.','DisplayName','dvx pre');
plot(real(analysis.pre.dvyKx) ./analysis.kyWavenums,'b-.','DisplayName','dvy pre');
plot(real(analysis.pre.dbyKx) ./analysis.kyWavenums,'c-.','DisplayName','dby pre');
plot(real(analysis.linear.omega(1:analysis.nModes(1)))./analysis.kyWavenums','k');
plot(0*numel(analysis.kyWavenums),'k-');

grid on;

xlabel('Mode number','fontsize',fsize);
ylabel('real(kx) / Ky','fontsize',fsize);

title('Fits of re(kx/ky) using multiple perturbed quantities','fontsize',fsize);


end
