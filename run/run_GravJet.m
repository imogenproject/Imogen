% Run a gravitational jet test

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run                 = JetInitializer();
run.mode.gravity    = true;
run.direction       = JetInitializer.X;
run.image.interval	= 3;
run.image.mass		= true;
run.image.mach		= true;
run.info            = 'Gravity jet test.';
run.notes           = '';

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
