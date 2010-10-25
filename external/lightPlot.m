function lightPlot(hFig)

    if ( nargin < 1 || isempty(hFig) ); hFig = get(0,'CurrentFigure'); end
    
    view(3); axis tight;
    camlight('headlight');
    material shiny;
    lighting gouraud;
    set(hFig,'Color','white','Renderer','OpenGL');
    set(hFig,'KeyPressFcn',@moveLightCB);
   
    hAxis = get(hFig,'CurrentAxes');
    hPlots = get(hAxis,'Children');
    for i=1:length(hPlots)
       if strcmpi(get(hPlots(i),'Type'),'surface')
         set(hPlots(i),'LineStyle','none','FaceColor','interp'); 
       end
    end
   
    %--- Context Menu ---%
    cmenu = figuresContextMenu(hFig, '3d');
    set(hAxis,'UIContextMenu',cmenu);
    for i=1:length(hPlots); set(hPlots(i),'UIContextMenu',cmenu); end
    
end

function moveLightCB(src, event)

    set(0,'CurrentFigure',src);
    
    switch (event.Key)
        case '='
            camlight('headlight');
        otherwise
            Colormap_KeyRelease_Callback(src, event);
    end

end