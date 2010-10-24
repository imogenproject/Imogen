function slicePlotter(data, pTag, slice)
% A visualization function to plot the intermediate slice results for 1D, 2D or 3D slice data.
%
%>> data	input data structure												struct
%>> pTag	tag of desired or existing plot										str
%>> slice	cutting indices for each grid direction								int   [sx sy sz]
%---------------------------------------------------------------------------------------------------
    
    %--- Initialization ---%
	if ( (nargin < 1 ) || isempty(data) )
		error('ImogenError:SlicePlotter','Missing data structure argument. Plotting aborted.');
	end
	
    if ( (nargin < 2) || isempty(pTag) ) %HANDLE: missing pTag arg
        for i=1:1000
            pTag = ['SlicePlot_' num2str(i)];
            if isempty(findobj('Tag',pTag)); break; end
        end
    end
    
    if ( (nargin < 3) || isempty(slice) ), slice = []; end
    

	%--- Determine and adjust array dimensions ---%
    dimension = ndims(data.mass);
    
    if (dimension < 3 && any(size(data.mass) == 1))
        dimension = 1;
        N = size(data.momX);
        if (N(1) == 1); data.momX = data.momX'; data.momY = data.momY'; data.momZ = data.momZ'; end
        N = size(data.magX);
        if (N(1) == 1); data.magX = data.magX'; data.magY = data.magY'; data.magZ = data.magZ'; end
        N = size(data.mass); 
        if (N(1) == 1); data.mass = data.mass'; data.ener = data.ener'; end
    end
    
    %--- Convert to run arrays ---%
    fields = {'X', 'Y', 'Z'};
    run  = fauxImogenManager(data.dGrid);
    
    mass = FluidArray(ENUM.SCALAR, ENUM.MASS, data.mass, run, []);
    ener = FluidArray(ENUM.SCALAR, ENUM.ENER,	data.ener, run, []);
	grav = GravityArray(ENUM.GRAV, run, []);
    grav.array = data.grav;
	
	mom	 = FluidArray.empty(3,0);
	mag	 = MagnetArray.empty(3,0);
	for i=1:3
		mom(i) = FluidArray(ENUM.VECTOR(i), ENUM.MOM, data.(['mom' fields{i}]), run, []);
		mag(i) = MagnetArray(ENUM.VECTOR(i), ENUM.MAG, data.(['mag' fields{i}]), run, []);
	end
    
	%--- Calculate Array Values ---%
	data.spen           = ener.array ./ mass.array;
	data.preT           = pressure('total', mass, mom, ener, mag, data.gamma);
	data.preG           = pressure('gas', mass, mom, ener, mag, data.gamma);
	data.preM           = pressure('magnet', mass, mom, ener, mag, data.gamma);
	data.delB           = magneticDivergence(run,mag);
	data.speed          = sqrt(getVelocitySquared(mass,mom));
	data.mach           = getMach(mass, mom, ener, mag, data.gamma);
	data.comp           = getCompression(run, mom, mass);
	data.vort           = getVorticity(run, mom, mass);
    data.magStrength    = getMagneticSquared(mag);
	data.beta           = data.preG ./ (0.5*max(data.magStrength, 1e-5));
    
    %--- Handle 3D Input ---%
    N = size(data.mass);
    guiInfo.activate = false;
    guiInfo.size = N;
    if (length(N) > 2 && isempty(data) )
        if isempty(slice) % Create default slice
            slice = ceil(N / 2);
            [minVal, index] = min(slice);
            slice(index) = 0;
        end
        hGUI = findobj('Tag',[pTag '_GUI']);
        if isempty(hGUI)
            guiInfo.activate = true;
            data3D = {'info', info;     'ver', ver;     'time', time; ...
                      'mass', mass;     'mom', mom;     'energy', energy;     'mag', mag; ...
                      'spen', spen;     'preT', preT;   'preG', preG;       'preM', preM; ...
                      'delB', delB;     'speed', speed;  'mach', mach;      'velComp', velComp; ...
                      'velVort', velVort;                'magStrength', magStrength; ...
                      'beta', beta ...
                      }; 
        end
    end
    
    %--- Slice the arrays if desired ---%
    if ~isempty(slice)
        sliceStr = local_sliceStrFromSlice(slice);
        mass    = squeeze( eval(['data.mass('      sliceStr]));
        energy  = squeeze( eval(['data.ener('    sliceStr]));
        
        mom     = squeeze( eval(['mom(:,'     sliceStr]));
        xMom    = squeeze( eval(['data.momX('      sliceStr]));
        yMom    = squeeze( eval(['data.momY('      sliceStr]));
        zMom    = squeeze( eval(['data.momZ('      sliceStr]));

        mag     = squeeze( eval(['mag(:,'      sliceStr]));
        xMag    = squeeze( eval(['data.magX('      sliceStr]));
        yMag    = squeeze( eval(['data.magY('      sliceStr]));
        zMag    = squeeze( eval(['data.magZ('      sliceStr]));

        spen    = squeeze( eval(['spen('      sliceStr]));
        preT    = squeeze( eval(['preT('      sliceStr]));
        preG    = squeeze( eval(['preG('      sliceStr]));
        preM    = squeeze( eval(['preM('      sliceStr]));
        delB    = squeeze( eval(['delB('      sliceStr]));
        speed   = squeeze( eval(['speed('     sliceStr]));
        mach    = squeeze( eval(['mach('      sliceStr]));
        velComp = squeeze( eval(['velComp('   sliceStr]));
        velVort = squeeze( eval(['velVort(:,' sliceStr]));
        magStrength = squeeze( eval(['magStrength(' sliceStr]));
        %beta    = squeeze( eval(['beta('      sliceStr]));
        
        N = size(mass);
    end
    
    %--- Determine Dimension of the Plots ---%
    if prod(double(N > 1)), dimension = true;
    else            dimension = false; 
    end
    
    %-----------------------------------------------------------------------------------------------
    % Plot data
    %----------
    
    %--- Find or create plot ---%
    hFig = findobj('Tag',pTag);
    if ~isempty(hFig); figure(hFig);
    else hFig = figure('Tag',pTag,'Name',pTag,'Color','white','KeyReleaseFcn',@Colormap_KeyRelease_Callback);
    end
    hold all;
    cmenu = figuresContextMenu(hFig,'2d');
    set(hFig,'UIContextMenu',cmenu);
    
    %--- Determine appropriate subplot arrangement ---%
    [rows, cols] = local_testForMagField(data.magX,data.magY,data.magZ);
    
    %--- Mass Plot ---%
    hAxis = subplot(rows,cols,1); hold all;
    hPlot = local_plot1or2D(data.mass,dimension,'mass', hAxis);

    %--- Total Energy Plot ---%
    hAxis = subplot(rows,cols,2); hold all;
    local_addHoldsArray(hAxis,'ETot','Total Energy',data.ener);
    local_addHoldsArray(hAxis,'Spen','Spec Energy',data.spen);
    hPlot = local_plotChoices(dimension,hAxis,hFig, true);

    %--- Pressure Plot ---%
    hAxis = subplot(rows,cols,3); hold all;
    local_addHoldsArray(hAxis,'preT','Total Press.',data.preT);
    local_addHoldsArray(hAxis,'preG','Gas Press.',data.preG);
    local_addHoldsArray(hAxis,'preM','Mag Press.',data.preM);
    hPlot = local_plotChoices(dimension,hAxis,hFig, true);

    %--- Gravity Plot ---%
    hAxis = subplot(rows,cols,4); hold all;
    local_addHoldsArray(hAxis,'potent','Gravity Pot.',zeros(size(data.mass)));
    hPlot = local_plotChoices(dimension,hAxis,hFig, true);
    
    %--- X Velocity Plot ---%
    hAxis = subplot(rows,cols,5); hold all;
    local_addHoldsArray(hAxis,'xvel','Vx', data.momX ./ data.mass);
    local_addHoldsArray(hAxis,'xvort','Vort X',squeeze(data.vort(1,:,:,:)));
    hPlot = local_plotChoices(dimension,hAxis,hFig, true);
    
    
    %--- Y Velocity Plot ---%
    hAxis = subplot(rows,cols,6); hold all;
    local_addHoldsArray(hAxis,'yvel','Vy', data.momY ./ data.mass);
    local_addHoldsArray(hAxis,'yvort','Vort Y',squeeze(data.vort(2,:,:,:)));
    hPlot = local_plotChoices(dimension,hAxis,hFig, true);
    
    %--- Z Velocity Plot ---%
    hAxis = subplot(rows,cols,7); hold all;
    local_addHoldsArray(hAxis,'zvel','Vz',data.momZ./ data.mass);
    local_addHoldsArray(hAxis,'zvort','Vort Z',squeeze(data.vort(3,:,:,:)));
    hPlot = local_plotChoices(dimension,hAxis,hFig, true);
    
    %--- Speed Plot ---%
    hAxis = subplot(rows,cols,8); hold all;
    local_addHoldsArray(hAxis,'speed','Speed',data.speed);
    local_addHoldsArray(hAxis,'mach','Mach',data.mach);
    local_addHoldsArray(hAxis,'comp','Compression',data.comp);
    hPlot = local_plotChoices(dimension,hAxis,hFig, true);
    
    %-----------------------------------------------------------------------------------------------
    % Plot Magnetic field if ~= 0
    %----------------------------
    if ( rows*cols >= 12 )
        %--- X Magnetic Field Plot ---%
        hAxis = subplot(rows,cols,9); hold all;
        hPlot = local_plot1or2D(data.magX,dimension,'X Mag.', hAxis);
        
        %--- Y Magnetic Field Plot ---%
        hAxis = subplot(rows,cols,10); hold all;
        hPlot = local_plot1or2D(data.magY,dimension,'Y Mag.', hAxis);

        %--- Z Magnetic Field Plot ---%
        hAxis = subplot(rows,cols,11); hold all;
        hPlot = local_plot1or2D(data.magZ,dimension,'Z Mag.', hAxis);
        
        %--- Magnetic Properties ---%
        hAxis = subplot(rows,cols,12); hold all;
        local_addHoldsArray(hAxis,'delB','Div.B',data.delB);
        local_addHoldsArray(hAxis,'strength','B Strength',data.magStrength);
        local_addHoldsArray(hAxis,'beta','Beta',data.beta);
        hPlot = local_plotChoices(dimension,hAxis,hFig, true);
    end
    
    %-----------------------------------------------------------------------------------------------
    %Label the graph with an overlay axis with text
    %----------------------------------------------
    loc = findstr(data.about, '\n'); 
    if isempty(loc); loc = length(info); else loc = loc(1) - 1; end
    try
        infoTime = time.time;
        infoIter = time.iteration;
    catch ME
        try 
            infoTime = time(2); 
            infoIter = time(1);
        catch ME
            infoTime = -1;
            infoIter = -1;
            warning('ImogenUtils:TimeData','Unable to determine time information for annotation.');
        end
    end
    txt = sprintf(['[Imogen v%s] [Iteration: %g at Time: %0.3g]: ' info(1:loc)],ver,infoIter,infoTime);
    hAno = annotation(hFig,'textbox',[0 0 1 0.04],'FontSize',8,'String',txt,'Interpreter','none');
    
    %--- Activate the GUI if needed ---%
    if (guiInfo.activate)
        hGUI = Slice3D_GUI(hFig, 'slicePlotter3x3', guiInfo.size, slice, data3D);
        set(hFig,'CloseRequestFcn',{@partneredCloseCB,hGUI});
        set(0,'CurrentFigure',hFig);
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cols = local_Set1DColors(currentPlotIndex)
    cols = zeros(10,3);
    cols(1,:) = [0 0 0];        %Black
    cols(2,:) = [1 0 0];        %Red
    cols(3,:) = [0 0 1];        %Blue
    cols(4,:) = [0.5 0 1];      %Purple
    cols(5,:) = [0 0.6 0];      %Green
    cols(6,:) = [1.0 0.5 0];    %Orange
    cols(7,:) = [0.6 0.6 0];    %Yellow
    cols(8,:) = [1.0 0.5 0.5];  %Pink
    cols(9,:) = [0.25 0.5 1];   %Aqua
    cols(10,:)= [0.5 0.75 0];   %Green-Yellow
    
    if ( nargin>0 && ~isempty(currentPlotIndex) ), cols = cols(currentPlotIndex,:); end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hPlot = local_plot1or2D(array,dimension,name, hAxis, currentPlotIndex)
    if ( nargin<5 || isempty(currentPlotIndex) ), currentPlotIndex = 1; end
    local_setAxisProps(hAxis,size(array), dimension);
    
    if (dimension)
        hPlot = imagesc(array');
        set(hPlot,'UserData',array');
    else    
        x = linspace(0,1,length(array));
        hPlot = plot(x,array,'LineWidth',2.0,'Color',local_Set1DColors(currentPlotIndex));
        set(hPlot,'UserData',array);
    end
    set(hPlot,'ButtonDownFcn',{@plotClick_Callback,name});
    ylabel(name);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function local_addHoldsArray(hAxis, field, name, array)
    holds = get(hAxis,'UserData');
    if isempty(holds); 
        holds.current = 0; 
        holds.list = {field}; 
    else
        N = length(holds.list);
        holds.list{N+1} = field;
    end
    minVal = min(min(min(min(array))));
    maxVal = max(max(max(max(array))));
    if (minVal == maxVal); maxVal = maxVal + 0.0001; end
    eval(['holds.' field '.array = array;']);
    eval(['holds.' field '.name = name;']);
    eval(['holds.' field '.min = minVal;']);
    eval(['holds.' field '.max = maxVal;']);
    
    set(hAxis,'UserData',holds);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hPlot = local_plotChoices(dimension,hAxis,hFig,ini)
    holds = get(hAxis,'UserData');
    N = length(holds.list);
    holds.current = holds.current + 1; if (holds.current > N); holds.current = 1; end
    
    set(hFig,'CurrentAxes',hAxis);
    curr = holds.list{holds.current};
    hPlot = local_plot1or2D(      eval(['holds.' curr '.array;']), ...
                            dimension,  eval(['holds.' curr '.name;']), hAxis, holds.current);
    set(hAxis,'CLim',[eval(['holds.' curr '.min;']) eval(['holds.' curr '.max;'])]);
    set(hAxis,'UserData',holds);
    if (ini); set(hAxis,'ButtonDownFcn',{@local_switchPlotCB,dimension,hFig}); end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function local_switchPlotCB(src,event,dimension,hFig)
    hs = get(src,'Children');
    for i=1:length(hs)
        type = get(hs(i),'Type');
        if ( strcmpi(type,'line') || strcmpi(type,'image') )
            delete(hs(i));
        end
    end
    hPlot = local_plotChoices(dimension,src,hFig,false);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function local_setAxisProps(hAxis, N, dimension)
    cols = local_Set1DColors(); %Customized color multiple 1D overlay plotting
    set(hAxis,'ColorOrder',cols,'XLim',[0 1],'FontSize',12);
    if (dimension) 
        M = N;
        [val,index] = min(N);
        if (length(M) == 2); M(3) = 1; end
        set(hAxis,'XLim',[1 M(1)],'YLim',[1 M(2)],'PlotBoxAspectRatio',M/N(index)); 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rows, cols] = local_testForMagField(xMag,yMag,zMag)
    rows = 3; cols = 4;
    magTest = [true true true];
    minVal = min(min(min(xMag))); maxVal = max(max(max(xMag))); if ( (minVal == 0) && (maxVal == 0) ) magTest(1) = false; end
    minVal = min(min(min(yMag))); maxVal = max(max(max(yMag))); if ( (minVal == 0) && (maxVal == 0) ) magTest(2) = false; end
    minVal = min(min(min(zMag))); maxVal = max(max(max(zMag))); if ( (minVal == 0) && (maxVal == 0) ) magTest(3) = false; end
    if ( ~magTest(1) && ~magTest(2) && ~magTest(3) ); rows = 3; cols = 3; end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sliceStr = local_sliceStrFromSlice(slice)
%- sliceStr = local_sliceStrFromSlice(slice) -------------------------------------------------------
% This routine will convert a slice array [sx sy sz] into a slice string (sx,sy,sz) for use in
% evaulating arrays at various slice points. (e.g. [0 12 5] -> :,12,5); ).
%---------------------------------------------------------------------------------------------------

    N = length(slice);
    sliceStr = '';
    for i=1:N % For each dimension in slice
        switch (slice(i))
            case 0 % Zero means all (i.e. 0 -> :)
                if isempty(sliceStr), sliceStr = ':';
                else                  sliceStr = [sliceStr ',:']; 
                end
            otherwise % Add number to slice
                if isempty(sliceStr), sliceStr = num2str(slice(i));
                else                  sliceStr = [sliceStr ',' num2str(slice(i))]; 
                end
        end
    end
    sliceStr = [sliceStr ');'];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function resArray = fix1DArray(inArray,isVector)
    N = size(inArray);
    resArray = inArray;
    if isVector
        
    else
        if (N(1) == 1 && length(N) < 3),     resArray = inArray'; end
    end

end
