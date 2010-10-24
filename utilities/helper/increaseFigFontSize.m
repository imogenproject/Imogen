function increaseFigFontSize(src)

    if ( (nargin < 1) || isempty(src) )
        src = get(0,'CurrentFigure');
    end

    hChildren = get(src, 'Children');
    
    for i=1:length(hChildren)
       type = get(hChildren(i), 'Type');
       disp(type);
       if ( strcmpi(type,'axes') || strcmpi(type,'text') )
            changeFontSize(hChildren(i));
            if strcmpi(type,'axes') 
                changeFontSize(get(hChildren(i), 'Title'));
                changeFontSize(get(hChildren(i), 'XLabel'));
                changeFontSize(get(hChildren(i), 'YLabel'));
                changeFontSize(get(hChildren(i), 'ZLabel'));
                increaseFigFontSize(hChildren(i));
            end
       end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function changeFontSize(src)
    fontSize = get(src, 'FontSize');
    fontSize = fontSize + 4;
    set(src,'FontSize', fontSize);
end
