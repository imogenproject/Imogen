function plotModePlane(dq, x, t, ky, kz)

fsize = 24;
subplot(1,2,1);

surf(t, x, log(abs(squeeze(dq(ky, kz, :, :)))),'linestyle','none');
xlabel('Distance from shock','fontsize',fsize);
ylabel('Time elapsed');
zlabel('ln |dq|');

subplot(1,2,2)
surf(t, x, unwrap(unwrap(squeeze(angle(dq(ky,kz,:,:))) ,1,2),1,1),'linestyle','none');
xlabel('Distance from shock','fontsize',fsize);
ylabel('Time elapsed');
zlabel('Mode''s phase');


end
