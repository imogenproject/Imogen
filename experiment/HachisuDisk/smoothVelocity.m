function mom = smoothVelocity(rho, momentum)

vel = momentum ./ rho;
vel(isnan(vel)) = 0;

vpr = vel;
vpr(vpr > 0) = 1;

for t = 1:7;
    lap = del2(vpr,1);
    lap(vpr == 0) = 0;
    vpr=vpr+lap;
end

%vpr(:,1)'
%exp((vpr(:,1)-1)/.5)'

mom = momentum .* exp((vpr-1)/.5);

end
