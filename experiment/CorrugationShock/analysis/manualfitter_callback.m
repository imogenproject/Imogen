function manualfitter_callback(src, eventdata, setterfcn, memfcn, analyzer, memory)

% memory:
% .ky
% .kz = the mode numbers we're working on.

% .kx      = the stored kx value guess
% .w       = the stored omega value guess
% .varfit  = which variable
% .typefit = which coefficient (kw and w, rel and imaginary)
% .qty     = 1 for post, 0 for pre.
% .dq      = change per key press in variables

% Keystrokes:
% [12345]   : [rho vx vy bx by] manual fit
% [zxcv]    : fit [wim wre kxim kxre]
% left-right: tilt line
% up-down   : move line up/down
% [asdw]    : change mode [-ky +ky -kz +kz]
% e         : toggle between pre and postshock quantities
% h: print help

memory = analyzer.manfit_state;

modechanged = 0;

shiftDown = any(strcmp(eventdata.Modifier,'shift'));
ctrlDown  = any(strcmp(eventdata.Modifier,'control'));
movefactor = 1 * (1+9*shiftDown)*(1+9*ctrlDown);

if strcmp(eventdata.Key,'1'); memory.varfit = 1; end
if strcmp(eventdata.Key,'2'); memory.varfit = 2; end
if strcmp(eventdata.Key,'3'); memory.varfit = 3; end
if strcmp(eventdata.Key,'4'); memory.varfit = 4; end
if strcmp(eventdata.Key,'5'); memory.varfit = 5; end

if strcmp(eventdata.Key,'z'); memory.typefit = 1; end
if strcmp(eventdata.Key,'x'); memory.typefit = 2; end
if strcmp(eventdata.Key,'c'); memory.typefit = 3; end
if strcmp(eventdata.Key,'v'); memory.typefit = 4; end

if strcmp(eventdata.Key,'a'); memory.ky = max(1, memory.ky - 1); modechanged = 1; end
if strcmp(eventdata.Key,'d'); memory.ky = min(analyzer.nModes(1), memory.ky+1); modechanged = 1; end
if strcmp(eventdata.Key,'s'); memory.kz = max(1, memory.kz - 1); modechanged = 1; end
if strcmp(eventdata.Key,'w'); memory.kz = min(analyzer.nModes(2), memory.kz+1); modechanged = 1; end

if strcmp(eventdata.Key,'e'); memory.qty = 1 - memory.qty; modechanged = 1; end

% L/R tilt the value
if strcmp(eventdata.Key,'leftarrow')
switch(memory.typefit) 
    case 1; memory.w(memory.varfit,1)  =  memory.w(memory.varfit,1) + 1i*memory.df(1)*movefactor;
    case 2; memory.w(memory.varfit,1)  =  memory.w(memory.varfit,1) +    memory.df(1)*movefactor;
    case 3; memory.kx(memory.varfit,1) = memory.kx(memory.varfit,1) + 1i*memory.df(1)*movefactor;
    case 4; memory.kx(memory.varfit,1) = memory.kx(memory.varfit,1) +    memory.df(1)*movefactor;
end
end

if strcmp(eventdata.Key,'rightarrow')
switch(memory.typefit) 
    case 1; memory.w(memory.varfit,1)  =  memory.w(memory.varfit,1) - 1i*memory.df(1)*movefactor;
    case 2; memory.w(memory.varfit,1)  =  memory.w(memory.varfit,1) -    memory.df(1)*movefactor;
    case 3; memory.kx(memory.varfit,1) = memory.kx(memory.varfit,1) - 1i*memory.df(1)*movefactor;
    case 4; memory.kx(memory.varfit,1) = memory.kx(memory.varfit,1) -    memory.df(1)*movefactor;
end 
end

% Shift the value up/down
if strcmp(eventdata.Key,'uparrow')
switch(memory.typefit) 
    case 1; memory.w(memory.varfit,2)  =  memory.w(memory.varfit,2) + 1i*memory.df(2)*movefactor;
    case 2; memory.w(memory.varfit,2)  =  memory.w(memory.varfit,2) +    memory.df(2)*movefactor;
    case 3; memory.kx(memory.varfit,2) = memory.kx(memory.varfit,2) + 1i*memory.df(2)*movefactor;
    case 4; memory.kx(memory.varfit,2) = memory.kx(memory.varfit,2) +    memory.df(2)*movefactor;
end 
end

if strcmp(eventdata.Key,'downarrow')
switch(memory.typefit)
    case 1; memory.w(memory.varfit,2)  =  memory.w(memory.varfit,2) - 1i*memory.df(2)*movefactor;
    case 2; memory.w(memory.varfit,2)  =  memory.w(memory.varfit,2) -    memory.df(2)*movefactor;
    case 3; memory.kx(memory.varfit,2) = memory.kx(memory.varfit,2) - 1i*memory.df(2)*movefactor;
    case 4; memory.kx(memory.varfit,2) = memory.kx(memory.varfit,2) -    memory.df(2)*movefactor;
end
end

if strcmp(eventdata.Key,'h')
fprintf('12345: Fit respectively for drho, dvx, dvy, dbx, dby\nzxcv : Fit im(w), re(w), im(kx), re(kx) respectively\nLeft/right to tilt line, up down to shift vertically.\nasdw to change mode numbers.\nq restores from oopsies\ne switches between post and preshock variables\n');
return;
end

if strcmp(eventdata.Key,'q') || (modechanged == 1) % quickfit
if memory.qty == 1
  memory.kx = [analyzer.post.drhoKx(memory.ky, memory.kz) analyzer.post.drhoK0(memory.ky, memory.kz); ...
               analyzer.post.dvxKx(memory.ky, memory.kz) analyzer.post.dvxK0(memory.ky, memory.kz); ...
               analyzer.post.dvyKx(memory.ky, memory.kz) analyzer.post.dvyK0(memory.ky, memory.kz); ...
               analyzer.post.dbxKx(memory.ky, memory.kz) analyzer.post.dbxK0(memory.ky, memory.kz); ...
               analyzer.post.dbyKx(memory.ky, memory.kz) analyzer.post.dbyK0(memory.ky, memory.kz) ];

  memory.w  = [analyzer.omega.fromdrho2(memory.ky, memory.kz) analyzer.omega.drho2_0(memory.ky, memory.kz); ...
               analyzer.omega.fromdvx2(memory.ky, memory.kz) analyzer.omega.dvx2_0(memory.ky, memory.kz); ...
               analyzer.omega.fromdvy2(memory.ky, memory.kz) analyzer.omega.dvy2_0(memory.ky, memory.kz); ...
               analyzer.omega.fromdbx2(memory.ky, memory.kz) analyzer.omega.dbx2_0(memory.ky, memory.kz); ...
               analyzer.omega.fromdby2(memory.ky, memory.kz) analyzer.omega.dby2_0(memory.ky, memory.kz)];


else
  memory.kx = [analyzer.pre.drhoKx(memory.ky, memory.kz) analyzer.pre.drhoK0(memory.ky, memory.kz); ...
               analyzer.pre.dvxKx(memory.ky, memory.kz) analyzer.pre.dvxK0(memory.ky, memory.kz); ...
               analyzer.pre.dvyKx(memory.ky, memory.kz) analyzer.pre.dvyK0(memory.ky, memory.kz); ...
               analyzer.pre.dbxKx(memory.ky, memory.kz) analyzer.pre.dbxK0(memory.ky, memory.kz); ...
               analyzer.pre.dbyKx(memory.ky, memory.kz) analyzer.pre.dbyK0(memory.ky, memory.kz) ];

  memory.w  = [analyzer.omega.fromdrho1(memory.ky, memory.kz) analyzer.omega.drho1_0(memory.ky, memory.kz); ...
               analyzer.omega.fromdvx1(memory.ky, memory.kz) analyzer.omega.dvx1_0(memory.ky, memory.kz); ...
               analyzer.omega.fromdvy1(memory.ky, memory.kz) analyzer.omega.dvy1_0(memory.ky, memory.kz); ...
               analyzer.omega.fromdbx1(memory.ky, memory.kz) analyzer.omega.dbx1_0(memory.ky, memory.kz); ...
               analyzer.omega.fromdby1(memory.ky, memory.kz) analyzer.omega.dby1_0(memory.ky, memory.kz)];
end
end % end autofit

memfcn(memory);
setterfcn(memory.ky, memory.kz, memory.w, memory.kx, memory.qty);
kx = memory.kx;
w  = memory.w;
rawline = [];
plotstyles = {'r','g','g--','b','b--'};

hold off;

ymin = 1e9; ymax = -1e9; % For axes clamp

for v = 1:5
    if memory.qty == 1
        if memory.typefit < 3
        switch(v)
            case 1; rawline = getWline(analyzer.post.drho, memory.ky, memory.kz, 1);
            case 2; rawline = getWline(analyzer.post.dvx,  memory.ky, memory.kz, 1);
            case 3; rawline = getWline(analyzer.post.dvy,  memory.ky, memory.kz, 1);
            case 4; rawline = getWline(analyzer.post.dbx,  memory.ky, memory.kz, 1);
            case 5; rawline = getWline(analyzer.post.dby,  memory.ky, memory.kz, 1);
        end
        else
        switch(v)
            case 1; rawline = getKxline(analyzer.post.drho, memory.ky, memory.kz, 1, analyzer.linearFrames((end-10):end));
            case 2; rawline = getKxline(analyzer.post.dvx,  memory.ky, memory.kz, 1, analyzer.linearFrames((end-10):end));
            case 3; rawline = getKxline(analyzer.post.dvy,  memory.ky, memory.kz, 1, analyzer.linearFrames((end-10):end));
            case 4; rawline = getKxline(analyzer.post.dbx,  memory.ky, memory.kz, 1, analyzer.linearFrames((end-10):end));
            case 5; rawline = getKxline(analyzer.post.dby,  memory.ky, memory.kz, 1, analyzer.linearFrames((end-10):end));
        end
        end

    else
        if memory.typefit < 3
        switch(v)
            case 1; rawline = getWline(analyzer.pre.drho, memory.ky, memory.kz, 0);
            case 2; rawline = getWline(analyzer.pre.dvx,  memory.ky, memory.kz, 0);
            case 3; rawline = getWline(analyzer.pre.dvy,  memory.ky, memory.kz, 0);
            case 4; rawline = getWline(analyzer.pre.dbx,  memory.ky, memory.kz, 0);
            case 5; rawline = getWline(analyzer.pre.dby,  memory.ky, memory.kz, 0);
        end
        else
        switch(v)
            case 1; rawline = getKxline(analyzer.pre.drho, memory.ky, memory.kz, 0, analyzer.linearFrames((end-10):end));
            case 2; rawline = getKxline(analyzer.pre.dvx,  memory.ky, memory.kz, 0, analyzer.linearFrames((end-10):end));
            case 3; rawline = getKxline(analyzer.pre.dvy,  memory.ky, memory.kz, 0, analyzer.linearFrames((end-10):end));
            case 4; rawline = getKxline(analyzer.pre.dbx,  memory.ky, memory.kz, 0, analyzer.linearFrames((end-10):end));
            case 5; rawline = getKxline(analyzer.pre.dby,  memory.ky, memory.kz, 0, analyzer.linearFrames((end-10):end));
        end
        end

    end

    xvals = {analyzer.pre.X, analyzer.post.X};
    plotqty = [];

    switch(memory.typefit)
        case 1;
            plotqty = mean(log(abs(rawline(:,analyzer.linearFrames))));
            plot(analyzer.frameTimes(analyzer.linearFrames), plotqty, plotstyles{v}, 'linewidth', 1+2*(v == memory.varfit) ); hold on;
        case 2;
            plotqty = mean(unwrap(angle(rawline(:,analyzer.linearFrames)),pi,2));
            plot(analyzer.frameTimes(analyzer.linearFrames), plotqty,  plotstyles{v}, 'linewidth', 1+2*(v == memory.varfit) ); hold on;
        case 3;
            plotqty = mean(log(abs(rawline)));
            plot(xvals{memory.qty+1}, plotqty, plotstyles{v}, 'linewidth', 1+2*(v == memory.varfit) ); hold on;
        case 4;
            plotqty = mean(unwrap(angle(rawline),pi,2));
            plot(xvals{memory.qty+1}, plotqty,  plotstyles{v}, 'linewidth', 1+2*(v == memory.varfit) ); hold on;
    end

    ymin = min(ymin, min(plotqty));
    ymax = max(ymax, max(plotqty));

    switch(memory.typefit)
        case 1; c = imag(memory.w(memory.varfit,:));
        case 2; c = real(memory.w(memory.varfit,:)); 
        case 3; c = imag(memory.kx(memory.varfit,:));
        case 4; c = real(memory.kx(memory.varfit,:));
    end

    plotvals={analyzer.frameTimes(analyzer.linearFrames), analyzer.frameTimes(analyzer.linearFrames), xvals{memory.qty+1},xvals{memory.qty+1} };

    plot(plotvals{memory.typefit}, plotvals{memory.typefit}*c(1) + c(2),'k-');
end

%ylim([1.1*min(ymin, c(2)) 1.1*max(ymax, c(2))]);

if memory.typefit < 3
    xlabel('Time','fontsize',16);
else
    xlabel('X','fontsize', 16);
end

if (memory.typefit == 1) || (memory.typefit == 3)
    ylabel('Amplitude','fontsize',16);
else
    ylabel('Phase','fontsize', 16);
end

fitstrings = {'drho','dvx','dvy','dbx','dby'};
fitvarstr = {'growth rate','frequency','damp rate','wavelength'};
sideid = {'preshock','postshock'};

sigma = [];
switch(memory.typefit)
    case 1; sigma = std(imag(memory.w(:,1)));
    case 2; sigma = std(real(memory.w(:,1)));
    case 3; sigma = std(imag(memory.kx(:,1)));
    case 4; sigma = std(real(memory.kx(:,1)));
end

title(sprintf('Mode [%i %i]: Fitting %s %s, %s\nCurrent sigma-hat: %f', memory.ky, memory.kz, fitstrings{memory.varfit}, fitvarstr{memory.typefit}, sideid{memory.qty+1}, sigma/analyzer.kyValues(memory.ky) ),'fontsize',16);


end % function

function linedat = getWline(dq, ky, kz, prepost)
if prepost == 1
    linedat = squeeze(dq(ky, kz, 2:20,:));
else
    linedat = squeeze(dq(ky, kz, (end-15):(end-1),:));
end
end

function linedat = getKxline(dq, ky, kz, prepost, tframes)
if prepost == 1
    linedat = squeeze(dq(ky, kz, :, tframes)).';
else
    linedat = squeeze(dq(ky, kz, :, tframes)).';
end

end

