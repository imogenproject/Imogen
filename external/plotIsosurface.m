function [facesAndVerts] = plotIsosurface(array,isoValue,saveObj,mode)

%=== PLOTISOSURFACE ================================================================================
%   This routine will create a 3D plot of iso (equal) valued surface for the supplied array using
%   the isoValue input as the surface value.
%===================================================================================================
%
%-CALLS---------------------------------------------------------------------------------------------
%   plotIsosurface(array,isoValue,mode)
%-INPUTs--------------------------------------------------------------------------------------------
% array     3D array of values to be plotted                                double  [Nx Ny Nz]
% isoValue  the value(s) on which to create the isosurface                  double  [n]
% saveName  name for an obj file for the isosurfaces                        bool    t/f
% mode      reserved input for later use                                    str     *
%---------------------------------------------------------------------------------------------------

    if (nargin < 3 || isempty(saveObj)), saveObj = false; end
    if (nargin < 4 || isempty(mode)), mode = 'single'; end
    mode = lower(mode);

    cols = jet(length(isoValue));
    
    hFig = figure();
    hold all;    
    switch (mode)
        case {'single','spread'}
            for i=1:length(isoValue)
                facesAndVerts = isosurface(array,isoValue(i));
                hPatch = patch(facesAndVerts); 
                set(hPatch,'FaceColor',cols(i,:),'EdgeColor','none','AmbientStrength',0);
                if (saveObj), saveObjMesh('', facesAndVerts); end
            end
        otherwise
    end
    
    view(3); axis tight;
    camlight('headlight');
    grid on;
    material shiny;
    lighting gouraud;
    set(hFig,'Color','white','Renderer','OpenGL');
    hAxis = get(hFig,'CurrentAxes');
    
    divider = 10;
    limLabels = {'XLim','YLim','ZLim'};
    for i=1:3
        lims = get(hAxis,limLabels{i});
        deltaLims = lims(2) - lims(1);
        set(hAxis,limLabels{i},[lims(1)-deltaLims/divider, lims(2)+deltaLims/divider]);
    end
    
    set(hAxis,'Box','on','AmbientLightColor',[0 0 0]);
    hold off;
    
end