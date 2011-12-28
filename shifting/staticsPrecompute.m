function [valstab coeffstab indextab] = staticsPrecompute(values, coeffs, indices, arrdims)
% Given a list of indices (x,y,z), the array dimensions, and associated values and coeffs,
% compute all possible permutations

nstats = size(values, 1);

valstab = zeros([nstats 6]);
coeffstab = zeros([nstats 6]);
indextab = zeros([nstats 6]);

% Permutation [x y z]
nx = arrdims(1); ny = arrdims(2); nz = arrdims(3);
linind = indices(:,1) + nx*(indices(:,2)-1) + nx*ny*(indices(:,3)-1);
[indextab(:,1) I] = sort(linind,1);
valstab(:,1)   = values(I);
coeffstab(:,1) = coeffs(I);

% Permutation [x z y]
nx = arrdims(1); ny = arrdims(3); nz = arrdims(2);
linind = indices(:,1) + nx*(indices(:,2)-1) + nx*ny*(indices(:,3)-1);
[indextab(:,2) I] = sort(linind,1);
valstab(:,2)   = values(I);
coeffstab(:,2) = coeffs(I);

% Permutation [y x z]
nx = arrdims(2); ny = arrdims(1); nz = arrdims(3);
linind = indices(:,2) + nx*(indices(:,1)-1) + nx*ny*(indices(:,3)-1);
[indextab(:,3) I] = sort(linind,1);
valstab(:,3)   = values(I);
coeffstab(:,3) = coeffs(I);


% Permutation [y z x]
nx = arrdims(2); ny = arrdims(3); nz = arrdims(1);
linind = indices(:,2) + nx*(indices(:,3)-1) + nx*ny*(indices(:,1)-1);
[indextab(:,4) I] = sort(linind,1);
valstab(:,4)   = values(I);
coeffstab(:,4) = coeffs(I);


% Permutation [z x y]
nx = arrdims(3); ny = arrdims(1); nz = arrdims(2);
linind = indices(:,3) + nx*(indices(:,1)-1) + nx*ny*(indices(:,2)-1);
[indextab(:,5) I] = sort(linind,1);
valstab(:,5)   = values(I);
coeffstab(:,5) = coeffs(I);


% Permutation [z y x]
nx = arrdims(3); ny = arrdims(2); nz = arrdims(1);
linind = indices(:,3) + nx*(indices(:,2)-1) + nx*ny*(indices(:,1)-1);
[indextab(:,6) I] = sort(linind,1);
valstab(:,6)   = values(I);
coeffstab(:,6) = coeffs(I);

% NOTE: this is correct for imogen, where Matlab indexes from 1
% NOTE: and we do not subtract 1 as is here in GPUimogen

end
