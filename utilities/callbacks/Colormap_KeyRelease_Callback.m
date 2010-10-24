function Colormap_KeyRelease_Callback(src, event)

    set(0,'CurrentFigure',src);
    
    switch (event.Key)
        case '1', cmap = jet(256);
        case '2', cmap = gray(256);
        case '3', cmap = bone(256);
        case '4', cmap = hot(256);
        case '5', cmap = hsv(256);
        case '6', cmap = copper(256);
        case '7', cmap = local_rgbvCol(256,'r');
        case '8', cmap = local_rgbvCol(256,'g');
        case '9', cmap = local_rgbvCol(256,'b');
        case '0', cmap = zeroCol(256);
            
        case 'q'
%             hs = get(src,'Children');
%             for i=1:length(hs)
%                 type = get(hs(i),'Type');
%                 if ( strcmpi(type,'axes') );
%                     hps = get(hs(i),'Children');
%                     for j=1:length(hps)
%                         hpType = get(hps(j),'Type');
%                         if ( strcmpi(hpType, 'image') )
%                             set(src,'CurrentAxes',hs(i));
%                             colorbar('location','southoutside');
%                             break
%                         end
%                     end
%                 end
%             end
            return
                    
        otherwise
            return
    end
    if any(strcmpi(event.Modifier,'shift'))
        cmap = flipdim(cmap,1);
    end
    
    set(src,'Colormap',cmap);
            
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cmap = local_rgbvCol(nCols,type)
    
    i = 1; j = 2; k = 3;
    switch (lower(type))
        case 'r'
            i = 1; j = 2; k = 3;
        case 'g'
            i = 3; j = 1; k = 2;
        case 'b'
            i = 2; j = 3; k = 1;
    end

    cmap = zeros(nCols, 3);
    delCol = floor([nCols/2 nCols/3]);
    
    cmap(:,k) = linspace(0,1,nCols);
    
    mic = delCol(1); mac = mic + delCol(1);
    cmap(mic:mac,i) = linspace(0,1,(mac-mic+1));
    
    mic = mac + 1; mac = nCols;
    cmap(mic:mac,j) = linspace(0,1,(mac-mic+1));
    
end