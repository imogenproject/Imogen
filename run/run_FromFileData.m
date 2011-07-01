%  Run a fluid jet test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
%        The number of X cells is the controlling variable in grid step size.
run                 = ExperimentFromFileInitializer([256 256 64]);
run.iterMax         = 10000;

run.fname = '/home/erik/newspheroid2.mat';
run.latheAboutZ = 0;
run.momFormat = 1; 

%run.rhoSphere;
%run.rhoBG;
%rho.enerSphere
run.mode.gravity = true;
run.gravity.constant    = 1;

run.bgDensityCoeff = 3e-5;

run.image.interval	= 1;
run.image.mass		= true;
run.image.speed		= true;
run.image.grav          = false;
run.image.ener          = false;

run.activeSlices.xy     = true;
run.activeSlices.yz     = true;
%run.activeSlices.xz     = false;
run.activeSlices.xyz 	= true;

run.ppSave.dim3 = 100*(100/10000); % every 20 steps
run.ppSave.dim2 = 10;

run.info            = 'First test of self gravitating toroid.';
run.notes           = '';

run.gravity.solver = 'biconj_gpu';
run.gravity.tolerance = 1e-10;
run.gravity.iterMax = 150;

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
