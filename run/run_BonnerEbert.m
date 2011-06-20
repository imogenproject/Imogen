%  Run a fluid jet test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run                 = BonnerEbertInitializer([96 96 96]);
run.iterMax         = 10;

%run.sphereK        = 1;  % Set pressure scaling constant
%run.sphereGAMMA    = 1;  % Set pressure scaling exponent
%run.rho0           = .1; % Set core density
%run.rmax           = 15; % Set radius of simulation domain

run.image.interval	= 1;
run.image.mass		= true;
run.image.speed		= true;
run.image.momX		= true;

run.activeSlices.xy     = true;
run.activeSlices.xyz 	= true;
run.ppSave.dim3 = 10;
run.ppSave.dim2 = 10;

run.info            = 'Bonner-Ebert sphere test.';
run.notes           = '';

run.gravity.solver = 'biconj';
run.gravity.tolerance = 1e-10;
run.gravity.iterMax   = 150;

%--- Run tests ---%
if (true)
    [mass, mom, ener, magnet, statics, ini] = run.getInitialConditions();
    IC.mass = mass;
    IC.mom = mom;
    IC.ener = ener;
    IC.magnet = magnet;
    IC.statics = statics;
    IC.ini = ini;
    icfile = [tempname '.mat'];

    save(icfile, 'IC');
    clear IC mass mom ener magnet statics ini run;
    imogen(icfile);
end

enderRun();
