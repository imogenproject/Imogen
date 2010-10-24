function [mass pos] = massQuantization(finestRho, dGrid, optPos)
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
nQuantLevels = max(ceil(log(grid)/log(2))) - 1;
mass         = cell(nQuantLevels, 1); % Contains an Nx by Ny by Nz array
pos          = cell(nQuantLevels, 1); % Contains a 3 by Nx by Ny by Nz array

%--- A note about which end the finest value is kept at ---%
%        At some point we end up counting both ways so it doesn't really matter.

%--- Get top-level positions ---%
%        If we get top-level positions, use them and multiply by the differential volume. If not,
%        go about computing positions.
if nargin == 3;
    pos{nQuantLevels} = optPos;
    mass{nQuantLevels} = finestRho .* dGrid{1} .* dGrid{2} .* dGrid{3}; 
else
    pos{nQuantLevels} = zeros([3 grid]);
    mass{nQuantLevels} = finestRho;
    % Generate initial centers for every grid location
    for n = 1:3
        if numel(dGrid{n}) > 1
            pos{nQuantLevels}(n,:,:,:) = cumsum(dGrid{n}, n) - .5*dGrid{n}(1,1,1);
            mass{nQuantLevels} = mass{nQuantLevels} .* dGrid{n}; % Compute mass = mass * dv 
        else
            mass{nQuantLevels} = mass{nQuantLevels} * dGrid{n};
            switch n
            case 1; b = (1:grid(1))'*dGrid{1}; pos{nQuantLevels}(1,:,:,:) = repmat(b, [1 grid(2) grid(3)]) - .5*dGrid{1};
            case 2; b = (1:grid(2))*dGrid{2}; pos{nQuantLevels}(2,:,:,:) = repmat(b, [grid(1) 1 grid(3)]) -.5*dGrid{2};
            case 3; b = rand([1 1 grid(3)]); b(1:grid(3)) = (1:grid(3))*dGrid{3}; pos{nQuantLevels}(3,:,:,:) = repmat(b, [grid(1) grid(2) 1]) - .5*dGrid{3};
            end
        end
    end
end

%--- Iteratively un-refine the grid to coarser and coarser levels ---%
oldgrid = grid;
for ct = (nQuantLevels-1):-1:1

    newgrid = ceil(oldgrid ./ 2);

    mass{ct} = zeros(newgrid);
    pos{ct} = zeros([3 newgrid]);

    for delx = 1:min(2, oldgrid(1)); for dely = 1:min(2, oldgrid(2)); for delz = 1:min(2, oldgrid(3))
        xbound = newgrid(1); if (delx == 2) && (mod(oldgrid(1), 2) == 1); xbound = max(xbound - 1, 1); end
        ybound = newgrid(2); if (dely == 2) && (mod(oldgrid(2), 2) == 1); ybound = max(ybound - 1, 1); end
        zbound = newgrid(3); if (delz == 2) && (mod(oldgrid(3), 2) == 1); zbound = max(zbound - 1, 1); end

        selx   = delx:2:oldgrid(1);
        sely   = dely:2:oldgrid(2);
        selz   = delz:2:oldgrid(3);

        mass{ct}(1:xbound, 1:ybound, 1:zbound) = mass{ct}(1:xbound, 1:ybound, 1:zbound) + ...
                                                 mass{ct+1}(selx,sely,selz);
    end; end; end


    for eta = 1:3;
    xps = zeros(newgrid);
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
