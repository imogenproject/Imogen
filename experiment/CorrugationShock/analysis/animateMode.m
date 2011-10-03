function animateMode(x, ky, kz, mode)

[xv yv] = ndgrid(x, 1:.1:6.28);

for t = 1:size(mode,4);
   dat = squeeze(mode(ky, kz, :, t));
   
   dat = dat*ones([1 size(yv,2)]);
   
   surf(log(abs(real(dat.*exp(1i*yv))')) );
   xlabel('X position')
   ylabel('Y position');
   title(t);
   pause(.2);
end


end