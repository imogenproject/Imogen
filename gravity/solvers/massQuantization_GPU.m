function [mass pos] = massQuantization_GPU(finestRho, dGrid, optPos)
% Computes a bi/quad/octree of mass and barycenters using the input finestRho and from integrating
% the given grid and dGrid values to find initial positions.
% 
%>> finestRho     Simulation-resolution array describing density              double
%>> dGrid          Grid step size, 3x1 cell array. Each cell may be either     double
%                  scalar (uniform spacing) or and array the same size as
%                  as finestRho (nonuniform spacing)
%>> optPos         Optional positions; Saving the returned finest-positio     double
%                  array and passing it back avoids re-integrating every step

grid = size(finestRho);
if numel(grid) == 2; grid(3) = 1; end

%--- Allocate cell arrays for quantization ---%
nQuantLevels = 0;
gridred = grid;
while (max(gridred - fix(gridred)) < 1e-8) && (min(gridred) >= 3)
    gridred = gridred/2;
    nQuantLevels = nQuantLevels+1;
end

nQuantlevels = nQuantLevels - 2;

mass         = cell(nQuantLevels, 1); % Contains an Nx by Ny by Nz array
pos          = cell(nQuantLevels, 1); % Contains a 3 by Nx by Ny by Nz array

pos0 = cell(3,1);

%--- A note about which end the finest value is kept at ---%
%        At some point we end up counting both ways so it doesn't really matter as long as we're consistent.

%--- Get top-level positions ---%
%        If we get top-level positions, use them and multiply by the differential volume. If not,
%        go about computing positions assigning <0,0,0> to the lower left corner
if nargin == 3;
    pos{nQuantLevels} = GPUdouble(optPos);
else
    for n = 1:3
        if numel(dGrid{n}) > 1
            pos0{n} = GPUdouble(cumsum(dGrid{n}, n) - .5*dGrid{n}(1,1,1));
        else
            switch n
            case 1; b = (1:grid(1))'*dGrid{1}; pos0{1} = GPUdouble(repmat(b, [1 grid(2) grid(3)]) - .5*dGrid{1});
            case 2; b = (1:grid(2))*dGrid{2}; pos0{2} = GPUdouble(repmat(b, [grid(1) 1 grid(3)]) -.5*dGrid{2});
            case 3; b = rand([1 1 grid(3)]); b(1:grid(3)) = (1:grid(3))*dGrid{3}; pos0{3} = GPUdouble(repmat(b, [grid(1) grid(2) 1]) - .5*dGrid{3});
            end
        end
    end

    pos{nQuantLevels} = zeros([3 size(pos0{1})], GPUdouble);
    pos{nQuantLevels}(1,:,:,:) = pos0{1};
    pos{nQuantLevels}(2,:,:,:) = pos0{2};
    pos{nQuantLevels}(3,:,:,:) = pos0{3};
end

mass{nQuantLevels} = GPUdouble(finestRho .* dGrid{1} .* dGrid{2} .* dGrid{3});

for ct = (nQuantLevels-1):-1:1
    % Store downsampled mass and position
    mass{ct} = interpolateGPUvar(mass{ct+1},-2);
    for x = 1:3; pos0{x} = interpolateGPUvar(pos0{x}.*mass{ct+1},-2); pos0{x} = pos0{x} ./ mass{ct}; end

    pos{ct} = zeros([3 size(pos0{1})], GPUdouble);
    pos{ct}(1,:,:,:) = pos0{1};
    pos{ct}(2,:,:,:) = pos0{2};
    pos{ct}(3,:,:,:) = pos0{3};   
end

return;

%--- Iteratively un-refine the grid to coarser and coarser levels ---%
oldgrid = grid;
for ct = (nQuantLevels-1):-1:1

    newgrid = ceil(oldgrid ./ 2);
    pos{ct} = zeros([3 newgrid],GPUdouble);

    %--- Loop over 3 dimensions [create one temp array instead of 3] ---%
    for eta = 1:3
        xps = zeros(newgrid,GPUdouble);
    
        %--- Count over the block of subcells ---%
        %        ... Which is between 2x2x2 and 1x1x1.
        for delx = 1:min(2, oldgrid(1)); for dely = 1:min(2, oldgrid(2)); for delz = 1:min(2, oldgrid(3))
            xbound = newgrid(1); if (delx == 2) && (mod(oldgrid(1), 2) == 1); xbound = max(xbound - 1, 1); end
            ybound = newgrid(2); if (dely == 2) && (mod(oldgrid(2), 2) == 1); ybound = max(ybound - 1, 1); end
            zbound = newgrid(3); if (delz == 2) && (mod(oldgrid(3), 2) == 1); zbound = max(zbound - 1, 1); end

            selx   = delx:2:oldgrid(1);
            sely   = dely:2:oldgrid(2);
            selz   = delz:2:oldgrid(3);

            xps(1:xbound,1:ybound,1:zbound)       = xps(1:xbound,1:ybound,1:zbound) + ...
                   mass{ct+1}(selx,sely,selz) .* squeeze(pos{ct+1}(eta,selx,sely,selz));

        end; end; end

        pos{ct}(eta,:,:,:) = xps ./ mass{ct};
    end

    oldgrid = newgrid;
end

end
