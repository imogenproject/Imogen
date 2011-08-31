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

rots = [];
fdir = [];
mdir = [];

switch(dim)
  case 2;
      switch( mod(run.time.iteration, 3) )
        case 1; rots = [0 2 0 0  0 2 0 0];
                 fdir = [1 2 0 0  2 1 0 0];
                 mdir = [2 1 0 0  3 3 0 0; ...
                         3 3 0 0  1 2 0 0];
        case 2; rots = [0 2 2 0  0 2 2 0];
                 fdir = [1 2 0 0  1 2 0 0];
                 mdir = [2 1 0 0  3 3 0 0; ...
                         3 3 0 0  2 1 0 0];
        case 0; rots = [2 2 0 0 2 2 0 0];
                 fdir = [2 1 0 0 2 1 0 0];
                 mdir = [1 2 0 0 3 3 0 0; ...
                         3 3 0 0 1 2 0 0 ];
      end
  case 3;
    switch( mod(run.time.iteration, 3) )
        case 1; rots = [0 2 3 0  0 2 3 0];
                 fdir = [1 2 3 0  3 2 1 0];
                 mdir = [2 1 1 0  2 3 3 0; ...
                         3 3 2 0  1 1 2 0];
        case 2; rots = [3 3 2 3  2 2 3 2];
                 fdir = [3 1 2 0  1 3 2 0];
                 mdir = [1 2 1 0  3 2 3 0; ...
                         2 3 3 0  2 1 1 0];
        case 0; rots = [2 3 2 2  3 2 3 3];
                 fdir = [2 3 1 0  2 1 3 0];
                 mdir = [1 1 2 0  3 3 2 0; ...
                         3 2 3 0  1 2 1 0];
    end
end

if order == 1; i0 = 1; else; i0 = 5; end

for i = i0:(i0+3)
%    if (rots(i) > 0) && (run.useGPU == true); xchgIndices(mass, mom, ener, mag, grav, rots(i)); end
    if run.fluid.ACTIVE && fdir(i) > 0;
        if run.useGPU == true; xchgIndices(mass, mom, ener, mag, grav, fdir(i)); end
	relaxingFluid(run, mass, mom, ener, mag, grav, fdir(i));
        if run.useGPU == true; xchgIndices(mass, mom, ener, mag, grav, fdir(i)); end
    end
    if run.magnet.ACTIVE && fdir(i) > 0; magnetFlux(run, mass, mom, mag, fdir(i) , mdir(:,i)); end
end

end

function xchgIndices(mass, mom, ener, mag, grav, toex)
l = [1 2 3];
l(1)=toex; l(toex)=1;

s = { mass, ener, mom(1), mom(2), mom(3) };

for i = 1:5
    s{i}.arrayIndexExchange(toex, 1);
    s{i}.store.arrayIndexExchange(toex, 0);
end

s = {mag(1).cellMag, mag(2).cellMag, mag(3).cellMag};
for i = 1:3
    s{i}.arrayIndexExchange(toex, 1);
%    s{i}.store.arrayIndexExchange(toex, 0);
end


end
