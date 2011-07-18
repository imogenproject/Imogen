function manualfitter_callback(src, eventdata, setterfcn, memfcn, analyzer, xaxis, yaxis, datain)

% memory:
% .kx = the stored kx value guess
% .w  = the stored omega value guess
% .df = the change if the arrows are pressed
% whofit = which curve we're trying to fit (toggle with space)

%eventdata
memory = analyzer.manfit_state;
%f strcmp(eventdata. 

if strcmp(eventdata.Key,'space')
    if memory.whofit == 1; memory.whofit = 0; else; memory.whofit = 1; end
end

% L/R tilt the w value (timewise)
if strcmp(eventdata.Key,'leftarrow')
    if memory.whofit == 1 % fit amplitude
        memory.kx = memory.kx - memory.df;
    else % fit phase
        memory.kx = memory.kx - 1i*memory.df;
    end
end

if strcmp(eventdata.Key,'rightarrow')
    if memory.whofit == 1 % fit amplitude
        memory.kx = memory.kx + memory.df;
    else % fit phase
        memory.kx = memory.kx + 1i*memory.df;
    end
end

% U/D tile the kx value (spacewise)
if strcmp(eventdata.Key,'uparrow')
    if memory.whofit == 1 % fit amplitude
        memory.w = memory.w - memory.df;
    else % fit phase
        memory.w = memory.w - 1i*memory.df;
    end
end

if strcmp(eventdata.Key,'downarrow')
    if memory.whofit == 1 % fit amplitude
        memory.w = memory.w + memory.df;
    else % fit phase
        memory.w = memory.w + 1i*memory.df;
    end
end

if strcmp(eventdata.Key,'w')
    if memory.whofit == 1
        memory.a0 = memory.a0 + memory.df;
    else
        memory.a0 = memory.a0 + 1i*memory.df;
    end
end

if strcmp(eventdata.Key,'s')
    if memory.whofit == 1
        memory.a0 = memory.a0 - memory.df;
    else
        memory.a0 = memory.a0 - 1i*memory.df;
    end
end

if strcmp(eventdata.Key,'a')
    memory.df = memory.df / 2;
    memfcn(memory); return;
end

if strcmp(eventdata.Key,'d')
    memory.df = memory.df * 2;
    memfcn(memory); return;
end

if strcmp(eventdata.Key,'h')
fprintf('Arrows: tilt plane\na/d: decrease/increase step size\nw/s: move plane down/up\nspace: work on amplitude/phase\n');
return;
end

if strcmp(eventdata.Key,'c')
setterfcn(memory.y, memory.z, memory.w, memory.kx, memory.qty)

end

memfcn(memory);

kx = memory.kx;
w  = memory.w;

hold off;
if memory.whofit == 1 % oscillatory part
    [u v] = ndgrid(xaxis, yaxis);
    kx = real(memory.kx);
    w = real(memory.w);

imgdat = unwrap(unwrap(angle(datain),2),1) - ( real(memory.a0) + kx*u + w*v);
    gain = 1;
    imgdat = imgdat*gain;
    imgdat(imgdat < -1) = -1;
    imgdat(imgdat > 1) = 1;

    surf(yaxis, xaxis, imgdat,'linestyle','none');
    hold on;
    contourf(yaxis, xaxis, imgdat , [0 0]);

    title(sprintf('Oscillatory component; a0=%g',real(memory.a0)));
else
    [u v] = ndgrid(xaxis, yaxis);
    kx = imag(memory.kx);
    w = imag(memory.w);

    imgdat =  log(abs(datain)) - ( imag(memory.a0) + kx*u + w*v);
    gain = 1;
    imgdat = imgdat*gain;
    imgdat(imgdat < -1) = -1;
    imgdat(imgdat > 1) = 1;

    surf(yaxis, xaxis, imgdat, 'linestyle','none');
    hold on;
    contourf(yaxis, xaxis, imgdat , [0 0]);
    title(sprintf('Amplitude component; a0=%g',imag(memory.a0)));
end

xlabel('time');
ylabel('distance');



end % function

