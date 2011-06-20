% Run a Sedov Taylor blast wave test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run         = SedovTaylorBlastWaveInitializer([65 65 65]);
run.iterMax = 100;

run.alias   = '';
run.info    = 'Sedov-Taylor blast wave test.';
run.notes   = 'Just a test...';



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
