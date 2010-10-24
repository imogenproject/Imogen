%   Run 3D Corrugation instability shock test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run          = CorrugationShockInitializer([300 48 48]);
run.profile  = true;
run.iterMax  = 2;
run.dataFile = 'ics-t10-ms10-ma0_5.mat';
run.info     = 'Corrugation instability test for parallel performance.';
run.notes    = 'Parallel performance test.';

%--- Run tests ---%
if (true) %Primary test
    [mass, mom, ener, magnet, statics, ini] = run.getInitialConditions();
    imogen(mass, mom, ener, magnet, ini, statics);
end

enderRun();
