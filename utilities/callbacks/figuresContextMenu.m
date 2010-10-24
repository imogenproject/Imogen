function cmenu = figuresContextMenu(src, mode)

    if ( (nargin < 1) || isempty(src) )
        src = get(0,'CurrentFigure');
    end

    if ( (nargin < 2) || isempty(mode) ) %HANDLE: missing mode arg
        mode = '2D';
    end
    mode = lower(mode);
    
    cmenu = uicontextmenu;
    
    if strcmpi(mode,'3d')
        item = uimenu(cmenu, 'Label', 'Update Light', 'Callback', ['camlight(''headlight'')']);
        
    end
    
    if any(strcmpi(mode,{'3d','2d'}))
        item = uimenu(cmenu, 'Label','Colormaps');
            colCB = sprintf('set(%s,''Colormap'',',num2str(src));
            c1 = uimenu(item, 'Label', 'Jet', 'Callback', [colCB 'jet(256))']);
            c2 = uimenu(item, 'Label', 'Gray', 'Callback', [colCB 'gray(256))']);
            c3 = uimenu(item, 'Label', 'Bone', 'Callback', [colCB 'bone(256))']);
            c4 = uimenu(item, 'Label', 'Hot', 'Callback', [colCB 'hot(256))']);
            c5 = uimenu(item, 'Label', 'HSV', 'Callback', [colCB 'hsv(256))']);
            c6 = uimenu(item, 'Label', 'Copper', 'Callback', [colCB 'copper(256))']);
            c7 = uimenu(item, 'Label', 'Zero', 'Callback', [colCB 'zeroCol(256))']);
    end
    
    pos = get(src,'Position');
    item = uimenu(cmenu, 'Label', 'Reset Size','Callback', ...
                   [sprintf('set(%s,''Position'',%s)',num2str(src),mat2str([pos(1),pos(2),560,420]))]);
              
    item = uimenu(cmenu, 'Label', 'Increase Font', 'Callback', [sprintf('increaseFigFontSize(%s)',num2str(src))]);
               
end