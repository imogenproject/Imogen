function plotClick_Callback(src, event, name)   
    array = get(src,'UserData');
    assignin('base','slicePlotterData',array);
    N = size(array); 
    dim2D = false; if ( (N(1) > 1) && (N(2) > 1) ); dim2D = true; end
    if (dim2D)
        answer = questdlg('Which plot type would you like?','Select Plot Type','Image','Surface','Cancel','Image');
        if strcmpi(answer,'Cancel'); return; end
        hFig = plot2DArray(array',100,answer);
    else
        hFig = figure('Color','white');
        hPlot = plot(linspace(0,1,length(array)),array,'LineWidth',2,'Color',get(src,'Color'));
        grid on;
    end
    hAxis = get(hFig,'CurrentAxes');
    hTitle = get(hAxis,'Title');
    set(hTitle,'String',name);
end