%   Run 3D Corrugation instability shock test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run         = CorrugationShockInitializer([1024 1024 1]);

run.iterMax     = 10000;
run.theta       = 0;
run.sonicMach   = 3;
run.alfvenMach  = .125;

run.ppSave.dim2 = .2;
run.ppSave.dim3 = 100;
run.seedAmplitude = 1e-5;

run.alias       = 'ms3_ma0pt125_gam53';
run.info        = sprintf('Corrugation instability test [Th=%g, Ms=%g, Ma=%g] with grid [%g, %g, %g]', ...
                          run.theta, run.sonicMach, run.alfvenMach, run.grid(1), run.grid(2), run.grid(3));
run.notes       = 'Corrugation instability test with maximal transverse resolution yet';

%--- Run tests ---%
if (true) %Primary test
    [mass, mom, ener, magnet, statics, ini] = run.getInitialConditions();
    imogen(mass, mom, ener, magnet, ini, statics);
end

enderRun();
