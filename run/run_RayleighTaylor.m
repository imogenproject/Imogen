%   Run 3D Corrugation instability shock test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run         = RayleighTaylorInitializer([128 384 1]);
run.info    = 'Rayleigh Taylor instability test';
run.notes   = '';

run.iterMax = 6000;

% run.rhoTop
% run.rhoBottom
% run.P0
% run.B0
% run.gravity.constant

%run.bcMode.gravity.y = 'const';

run.pertAmplitude = .02;
%run.Kx = 4;
run.randomPert = 1;
% run.Kz

run.image.interval = 5;
run.image.mass = true;
run.image.ener = true;
run.image.speed = true;
run.image.pGas = true;

run.ppSave.dim2 = 1;

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
