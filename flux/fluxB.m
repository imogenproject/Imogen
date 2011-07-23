function flux(run, mass, mom, ener, mag, grav, order)
% This function manages the fluxing routines for the split code by managing the appropriate fluxing
% order to average out any biasing caused by the Strang splitting.
%
%>< run run manager object ImogenManager
%>< mass mass density FluidArray
%>< mom momentum density FluidArray(3)
%>< ener energy density FluidArray
%>< mag magnetic field MagnetArray(3)
%>< grav gravitational potential GravityArray
%>> order direction of flux sweep (1 forward/-1 backward) int +/-1

dim = 3;
if (mass.gridSize(3) == 1); dim = 2; end
if (mass.gridSize(2) == 1) && (mass.gridSize(3) == 1); dim = 1; end

if dim == 1
    if run.fluid.ACTIVE; relaxingFluid(run, mass, mom, ener, mag, grav, 1); end
    %if run.magnet.ACTIVE; magnetFlux(run, mass, mom, mag, directVec(n), magneticIndices(n,:)); end
    return;
end

switch(dim)
  case 2;
      switch( mod(run.time.iteration, 3) )
        case 1; rots = [0 2 0 0  0 2 0 0];
                 fdir = [1 2 0 0  2 1 0 0];
                 mdir = [ ];
        case 2; rots = [0 2 2 0  0 2 2 0];
                 fdir = [1 2 0 0  1 2 0 0];
                 mdir = [ ];
        case 3; rots = [2 2 0 0 2 2 0 0];
                 fdir = [2 1 0 0 2 1 0 0];
                 mdir = [ ];
      end
  case 3;
    switch( mod(run.time.iteration, 3) )
        case 1; rots = [0 2 3 0  0 2 3 0];
                 fdir = [1 2 3 0  3 2 1 0];
                 mdir = [ ];
        case 2; rots = [3 3 2 3  2 2 3 2];
                 fdir = [3 1 2 0  1 3 2 0];
                 mdir = [ ];
        case 3; rots = [2 3 2 2  3 2 3 3];
                 fdir = [2 3 1 0  2 1 3 0];
                 mdir = [ ];
    end
end

if order == 1; i0 = 1; else; i0 = 5; end

for i = i0:(i0+3)
    if (rots(i) > 0) && (run.useGPU == true); xchgIndices(mass, mom, ener, mag, grav, rots(i)); end
    if run.fluid.ACTIVE && fdir(i) > 0; relaxingFluid(run, mass, mom, ener, mag, grav, fdir(i)); end
    %if run.magnet.ACTIVE; magnetFlux(run, mass, mom, mag, directVec(n), magneticIndices(n,:)); end
end

%    case 1
%     1 2 3 3 2 1
%     2 1 1 2 3 3
%     3 3 2 1 1 2
%    case 2
%     3 1 2 1 3 2
%     1 2 1 3 2 3
%     2 3 3 2 1 1
%    case 3
%     2 3 1 2 1 3
%     1 1 2 3 3 2
%     3 2 3 1 2 1

end

function xchgIndices(mass, mom, ener, mag, grav, toex)

s = { mass, ener, mom(1), mom(2), mom(3), mag(1), mag(2), mag(3) };

l = [1 2 3];
l(1)=toex; l(toex)=1;

for i = 1:8
    s{i}.arrayIndexExchange(toex, 1);
    if i < 6; s{i}.store.arrayIndexExchange(toex, 0); end
end

end
