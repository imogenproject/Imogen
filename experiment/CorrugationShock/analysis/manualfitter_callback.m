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
end

if strcmp(eventdata.Key,'d')
    memory.df = memory.df * 2;
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
    [u v] = ndgrid(xaxis, yaxis(analyzer.linearFrames));
    kx = real(memory.kx);
    w = real(memory.w);

imgdat = unwrap(unwrap(angle(datain),2,2),2,1);
    gain = 1;
    imgdat = imgdat*gain;

    subplot(1,2,1);
    [az,el] = view();
    hold off;
    surf(yaxis(analyzer.linearFrames), xaxis, imgdat(:,analyzer.linearFrames),'linestyle','none');
    hold on;
    surf(yaxis(analyzer.linearFrames), xaxis, real(memory.a0) + kx*u + w*v,'linestyle','none');
    view([az, el]);

    %contourf(yaxis, xaxis, imgdat , [0 0]);
    title(sprintf('Oscillatory component; a0=%g, df=%g',real(memory.a0), memory.df));
        
    subplot(1,2,2);
    contourf(yaxis(analyzer.linearFrames)', xaxis, imgdat(:, analyzer.linearFrames)-(real(memory.a0) + kx*u + w*v),[-.5 -.25 0 .25 .5]);
else
    [u v] = ndgrid(xaxis, yaxis(analyzer.linearFrames));
    kx = imag(memory.kx);
    w = imag(memory.w);

    imgdat =  log(abs(datain)) ;
    gain = 1;
    imgdat = imgdat*gain;
    %imgdat(imgdat < -1) = -1;
    %imgdat(imgdat > 1) = 1;

    hold off; 
    
    subplot(1,2,1);
    [az,el] = view();
    hold off;
    surf(yaxis(analyzer.linearFrames)', xaxis, imgdat(:, analyzer.linearFrames), 'linestyle','none');
    hold on;
    surf(yaxis(analyzer.linearFrames)', xaxis, imag(memory.a0) + kx*u + w*v,'linestyle','none');
    title(sprintf('Amplitude component; a0=%g, df=%g',imag(memory.a0), memory.df));
    view([az el]);    

    subplot(1,2,2);
    %surf(yaxis(analyzer.linearFrames)', xaxis, imgdat(:, analyzer.linearFrames)-(imag(memory.a0) + kx*u + w*v),'linestyle','none');
    contour3(yaxis(analyzer.linearFrames)', xaxis, imgdat(:, analyzer.linearFrames)-(imag(memory.a0) + kx*u + w*v));
end

xlabel('time');
ylabel('distance');



end % function


