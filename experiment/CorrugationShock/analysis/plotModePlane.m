function plotModePlane(dq, x, t, ky, kz)

fsize = 24;
subplot(1,2,1);

surf(t, x, log(abs(squeeze(dq(ky, kz, :, :)))),'linestyle','none');
ylabel('Distance from shock','fontsize',fsize);
xlabel('Time elapsed','fontsize',fsize);
zlabel('ln |dq|','fontsize',fsize);

subplot(1,2,2)
surf(t, x, unwrap(unwrap(squeeze(angle(dq(ky,kz,:,:))) ,1,2),1,1),'linestyle','none');
ylabel('Distance from shock','fontsize',fsize);
xlabel('Time elapsed','fontsize',fsize);
zlabel('Mode''s phase','fontsize',fsize);


end
