function [mass ener] = generateTestSystem(grid, rhosphere, rhobg, G, K0)

%--- Create some constants ---%
[X Y Z] = ndgrid(1:grid(1), 1:grid(2), 1:grid(3));
ballRad = .375;

%--- Create sphere on left-hand side ---%
R1 = 2*sqrt((X-grid(1)/4 + .5).^2 + (Y - grid(2)/2 + .5).^2 + (Z - grid(3)/2 + .5).^2) / grid(1);

pcRad = 0:1/grid(2):3;
  pcMass(pcRad < ballRad) = rhosphere;
  pcMass(pcRad >= ballRad)= rhobg;

  % The factor of 1.5: I don't think my setting of gamma to 2 from 5/3 is getting through
  pcEner(pcRad < ballRad) = rhosphere * (K0 + (2*pi*G*rhosphere/3)*(ballRad^2 - pcRad(pcRad < ballRad).^2));
  pcEner(pcRad > ballRad) = K0 * rhosphere; % Maintain equal pressure outside

mass = pchip(pcRad, pcMass, R1);
ener = pchip(pcRad, pcEner, R1);

%mass(R1 <  ballRad) = rhosphere;
%mass(R1 >= ballRad) = rhobg;

%ener(R1 >= ballRad) = K0 * rhosphere;

%mass = reshape(mass, grid);
%ener = reshape(ener, grid);

%--- Create mirror image on right hand side ---%
%        We may someday make a non-spherical object so mirror properly.
	mass(grid(1)/2 + 1:grid(1),:,:) = mass(grid(1)/2:-1:1,:,:);
	ener(grid(1)/2 + 1:grid(1),:,:) = ener(grid(1)/2:-1:1,:,:);

%--- Some work to avoid initial transients ---%
%mass = simpleBlur(mass, .6, 1);
%ener = simpleBlur(ener, .6, 1);


end
