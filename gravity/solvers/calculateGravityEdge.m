function massGrav = calculateGravityEdge(mass, dGrid, mirrorZ, bcsrc)
% Generates the right hand side for an iterative solver, consisting of both mass density (scaled by
% 4 pi G) and subtracting off potential at the boundaries (with scalings determined by the
% Laplacian stencil being used).
%
%>> mass      Mass density object                                              struct
%>> dgrid     3x1 cell array describing grid steps                             cell array
%<< massGrav  Rho + the boundary conditions, converted into an Nx1 vector      double

%--- Step 1: Calculate potential at boundaries ---%
grid = mass.gridSize;
N = grid+2;
coarseningConstant = sqrt(dGrid{1}^2 + dGrid{2}^2 + dGrid{3}^2)/sqrt(3);

[rhos poss] = massQuantization(mass.array, dGrid);

%--- Select boundary condition source ---%
if strcmp(bcsrc, ENUM.GRAV_BCSOURCE_FULL);   bcfinder = @mg_bc_matlab; end
if strcmp(bcsrc, ENUM.GRAV_BCSOURCE_INTERP); bcfinder = @interpolateBvec; end;

%--- Compute potential across 6 outer faces ---%
%        bvec = [-x -y -z +x +y +z nx-1 ny-1 nz-1]
%        Passing a simple array of cells to examine is better and could be done with a single call.

Xmin = -.5 * dGrid{1};
Xmax = (grid(1) + .5)*dGrid{1};
Ymin = -.5 * dGrid{2};
Ymax = (grid(2) + .5)*dGrid{2};
Zmin = -.5 * dGrid{3};
Zmax = (grid(3) + .5)*dGrid{3};
Nx = grid(1)+1;
Ny = grid(2)+1;
Nz = grid(3)+1;

if mirrorZ == false
    bvectorSet = cell([6 1]);
else
    bvectorSet = cell([12 1]);
end

%--- Faces of constant X ---%
bvectorSet{1} = [Xmin Ymin Zmin Xmin Ymax Zmax 0 Ny Nz];
bvectorSet{2} = [Xmax Ymin Zmin Xmax Ymax Zmax 0 Ny Nz];
if mirrorZ == true
    bvectorSet{7} = [Xmin Ymin -Zmax Xmin Ymax -Zmin 0 Ny Nz];
    bvectorSet{8} = [Xmax Ymin -Zmax Xmax Ymax -Zmin 0 Ny Nz];
end

%--- Faces of constant Y
bvectorSet{3} = [Xmin Ymin Zmin Xmax Ymin Zmax Nx 0 Nz];
bvectorSet{4} = [Xmin Ymax Zmin Xmax Ymax Zmax Nx 0 Nz];
if mirrorZ == true
    bvectorSet{9}  = [Xmin Ymin -Zmax Xmax Ymin -Zmin Nx 0 Nz];
    bvectorSet{10} = [Xmin Ymax -Zmax Xmax Ymax -Zmin Nx 0 Nz];
end

bvectorSet{5} = [Xmin Ymin Zmin Xmax Ymax Zmin Nx Ny 0];
bvectorSet{6} = [Xmin Ymin Zmax Xmax Ymax Zmax Nx Ny 0];
if mirrorZ == true
    bvectorSet{11} = [Xmin Ymin -Zmin Xmax Ymax -Zmin Nx Ny 0];
    bvectorSet{12} = [Xmin Ymin -Zmax Xmax Ymax -Zmax Nx Ny 0];
end

phiSet = bcfinder(rhos, poss, bvectorSet, coarseningConstant, dGrid{1});

if mirrorZ == false;
    lowerYZ = squeeze(phiSet{1});
    upperYZ = squeeze(phiSet{2});

    lowerXZ = squeeze(phiSet{3});
    upperXZ = squeeze(phiSet{4});

    lowerXY = squeeze(phiSet{5});
    upperXY = squeeze(phiSet{6});
else
    lowerYZ = squeeze(phiSet{1}+phiSet{7}(1,:,end:-1:1) );
    upperYZ = squeeze(phiSet{2}+phiSet{8}(1,:,end:-1:1) );

    lowerXZ = squeeze(phiSet{3}+phiSet{9}(:,1,end:-1:1) );
    upperXZ = squeeze(phiSet{4}+phiSet{10}(:,1,end:-1:1) );

    lowerXY = squeeze(phiSet{5}+phiSet{11});
    upperXY = squeeze(phiSet{6}+phiSet{12});
end

lowerXZ(1,:) = 0; lowerXZ(N(1),:) = 0;
upperXZ(1,:) = 0; upperXZ(N(1),:) = 0;

lowerXY(1,:) = 0; lowerXY(N(1),:) = 0;
lowerXY(:,1) = 0; lowerXY(:,N(2)) = 0;

upperXY(1,:) = 0; upperXY(N(1),:) = 0;
upperXY(:,1) = 0; upperXY(:, N(2)) = 0;

%--- 2. Apply to truncated parts of Laplacian stencil ---%
Nm = N - 1;
Nmm = Nm - 1;
massGrav = zeros(grid);

%--- Various gravity stencils; Only HOC6 is likely to be used ---%
%        For HOC6
%shifts = [0 0; 1 0; -1 0; 0 1; 0 -1; 1 1; 1 -1; -1 1; -1 -1];
%preFactors = [14, 3, 3, 3, 3, 1, 1, 1, 1];

%        For HOC4
shifts = [0 0; 1 0; -1 0; 0 1; 0 -1; 1 1; 1 -1; -1 1; -1 -1];
preFactors = [2, 1, 1, 1, 1, 0, 0, 0, 0];

%        Standard 7-point
%shifts = [0 0; 1 0; -1 0; 0 1; 0 -1; 1 1; 1 -1; -1 1; -1 -1];
%preFactors = [1, 0, 0, 0, 0, 0, 0, 0, 0];

%tic;
for i=1:9
    bcVals = circshift(lowerYZ,shifts(i,:)); 
    massGrav(1,:,:) = squeeze(massGrav(1,:,:)) + preFactors(i) * bcVals(2:Nm(2),2:Nm(3));
    bcVals = circshift(upperYZ,shifts(i,:)); 
    massGrav(Nmm(1),:,:) = squeeze(massGrav(Nmm(1),:,:)) + preFactors(i) * bcVals(2:Nm(2),2:Nm(3));
    
    bcVals = circshift(lowerXZ,shifts(i,:)); 
    massGrav(:,1,:) = squeeze(massGrav(:,1,:)) + preFactors(i) * bcVals(2:Nm(1),2:Nm(3));
    bcVals = circshift(upperXZ,shifts(i,:)); 
    massGrav(:,Nmm(2),:) = squeeze(massGrav(:,Nmm(2),:)) + preFactors(i) * bcVals(2:Nm(1),2:Nm(3));
    
    bcVals = circshift(lowerXY,shifts(i,:)); 
    massGrav(:,:,1) = squeeze(massGrav(:,:,1)) + preFactors(i) * bcVals(2:Nm(1),2:Nm(2));
    bcVals = circshift(upperXY,shifts(i,:)); 
    massGrav(:,:,Nmm(3)) = squeeze(massGrav(:,:,Nmm(3))) + preFactors(i) * bcVals(2:Nm(1),2:Nm(2));
end

rhsForce = .5*mass.array + (circ_shift(mass.array,1,1) + circ_shift(mass.array,1,-1) ...
                          +  circ_shift(mass.array,2,1) + circ_shift(mass.array,2,-1) ...
                          +  circ_shift(mass.array,3,1) + circ_shift(mass.array,3,-1))/12;

rhsForce(1,:,:) = rhsForce(1,:,:) - mass.array(end,:,:)/12;
rhsForce(end,:,:) = rhsForce(end,:,:) - mass.array(1,:,:)/12;

rhsForce(:,1,:) = rhsForce(:,1,:) - mass.array(:,end,:)/12;
rhsForce(:,end,:) = rhsForce(:,end,:) - mass.array(:,1,:)/12;

rhsForce(:,:,1) = rhsForce(:,:,1) - mass.array(:,:,end)/12;
rhsForce(:,:,end) = rhsForce(:,:,end) - mass.array(:,:,1)/12;

%--- 3. Add mass density and reshape into vector to complete right hand side ---%
massGrav = reshape(-massGrav/(6*dGrid{1}^2) + 4*pi*rhsForce, [prod(grid) 1]);

end

