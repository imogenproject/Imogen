function omegaAnalysis(analysis, fsize)

if nargin == 1; fsize = 16; end

hold off; 

subplot(1,2,1);

plot(-imag(analysis.post.drhoKx)./analysis.kyValues,'r','DisplayName','drho post');
hold on;
plot(-imag(analysis.post.dvxKx) ./analysis.kyValues,'g-','DisplayName','dvx post');
plot(-imag(analysis.post.dvyKx) ./analysis.kyValues,'b-','DisplayName','dvy post');
plot(-imag(analysis.post.dbxKx) ./analysis.kyValues,'m-','DisplayName','dbx post');
plot(-imag(analysis.post.dbyKx) ./analysis.kyValues,'c-','DisplayName','dby post');

plot(-imag(analysis.pre.drhoKx)./analysis.kyValues,'r-.','DisplayName','drho pre');
plot(-imag(analysis.pre.dvxKx) ./analysis.kyValues,'g-.','DisplayName','dvx pre');
plot(-imag(analysis.pre.dvyKx) ./analysis.kyValues,'b-.','DisplayName','dvy pre');
plot(-imag(analysis.pre.dbxKx) ./analysis.kyValues,'m-.','DisplayName','dbx pre');
plot(-imag(analysis.pre.dbyKx) ./analysis.kyValues,'c-.','DisplayName','dby pre');

plot(1:analysis.nModes(1), zeros(size(analysis.kyValues)),'k-','DisplayName','');

xlabel('Mode number','fontsize',fsize);
ylabel('imag(kx) / Ky','fontsize',fsize);

title('Fits of im(kx/ky) using multiple perturbed quantities','fontsize',fsize);

grid on

subplot(1,2,2);

plot(real(analysis.post.drhoKx)./analysis.kyValues,'r','DisplayName','drho post');
hold on;
plot(real(analysis.post.dvxKx) ./analysis.kyValues,'g-','DisplayName','dvx post');
plot(real(analysis.post.dvyKx) ./analysis.kyValues,'b-','DisplayName','dvy post');
plot(real(analysis.post.dbxKx) ./analysis.kyValues,'m-','DisplayName','dbx post');
plot(real(analysis.post.dbyKx) ./analysis.kyValues,'c-','DisplayName','dby post');

plot(real(analysis.pre.drhoKx)./analysis.kyValues,'r-.','DisplayName','drho pre');
plot(real(analysis.pre.dvxKx) ./analysis.kyValues,'g-.','DisplayName','dvx pre');
plot(real(analysis.pre.dvyKx) ./analysis.kyValues,'b-.','DisplayName','dvy pre');
plot(real(analysis.pre.dbxKx) ./analysis.kyValues,'m-.','DisplayName','dbx pre');
plot(real(analysis.pre.dbyKx) ./analysis.kyValues,'c-.','DisplayName','dby pre');

plot(1:analysis.nModes(1), zeros(size(analysis.kyValues)),'k-','DisplayName','');
grid on;

xlabel('Mode number','fontsize',fsize);
ylabel('real(kx) / Ky','fontsize',fsize);

title('Fits of re(kx/ky) using multiple perturbed quantities','fontsize',fsize);


end
