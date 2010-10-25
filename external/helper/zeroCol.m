function colors = zeroCol(nCols)

    if ( (nargin < 1) || isempty(nCols) ) %HANDLE: missing nColors arg
        nCols = 256;
    end
    
    hfm = floor(nCols/2); hfp = ceil(nCols/2); 
    if (hfm == hfp); hfp = hfm + 1; end
    colors = zeros(nCols,3);
    neg = linspace(1,0,hfm)';
    pos = linspace(0,1,hfm)';
    colors(1:hfm,1) = neg;
    for j=1:3; colors(hfp:nCols,j) = pos; end

end