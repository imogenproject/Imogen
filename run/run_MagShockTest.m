% Runs a magnetic shock tube problem.

%--- Initialize Imogen directory ---%
starterRun();


%--- Initialize test ---%
run               = MagneticShockTubeInitializer([1024 3 3]);

run.alias         = '';
run.info          = 'Magnetic shock tube test';
run.notes         = '';

%--- Run tests ---%
if (true) % Primary test
    run.alias  = [run.alias, '_BXBY'];
    run.xField = true;
    run.yField = true;
    run.zField = false;
    run.direction = MagneticShockTubeInitializer.X;
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

if (false) % Secondary test
    run.alias  = [run.alias, '_BXBZ'];
    run.xField = true;
    run.yField = false;
    run.zField = true;
    run.direction = MagneticShockTubeInitializer.X;
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
