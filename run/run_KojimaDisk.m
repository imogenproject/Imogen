%   Run a Kojima disk model.

%--- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
grid                = [384 384 1];
run                 = KojimaDiskInitializer(grid);
run.iterMax         = 1000;
run.save            = true;
run.edgePadding     = 0.5;
run.pointRadius     = 0.35;
run.radiusRatio     = 0.8;

run.image.interval  = 3;
run.image.speed     = true;
run.image.mass      = true;

run.specSaves.dim2  = 1:10;
run.ppSave.dim2     = 5;

run.bcMode.x        = ENUM.BCMODE_FADE;
run.bcMode.y        = ENUM.BCMODE_FADE;

run.addFade(ceil(grid/2), 16, ENUM.POINT_FADE , true, {ENUM.MOM});

run.info        = 'Kojima disk simulation';
run.notes       = '';

%--- Run tests ---%
if (true) %Primary test
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
