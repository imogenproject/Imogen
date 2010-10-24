function result = wormhole_shift(array, DIRECT, n, obj)
% Shifts an array circularly, so that the cells lost at one end are added as the new cells gained.
%
%>> array   array to be shifted                                             double((3),Nx,Ny,Nz)
%>> DIRECT  direction/dimension along which to shift                        int
%>> n       number of cells to shift                                        int
%>< obj     object owning the array                                         ImogenArray...
%<< result  shifted array                                                   double((3),Nx,Ny,Nz)


    %--- Initialization ---%
    N         = size(array);
    NSize     = length(N);
    index     = cell(1, NSize);
%    fadeSize  = obj.bcInfinity;
    for i=1:NSize,  index{i} = 1:N(i); end
    
    amt = abs(n);
    dir = sign(n);

    %--- Wormhole the - edge if we're shifting positively ---%
if dir == -1;
    wormIndex = index;
    wormIndex{DIRECT}    = 1:amt;

    previousDir = [3 1 2];
    nextDir     = [2 3 1];

    wormIndex{previousDir(DIRECT)} = index{previousDir(DIRECT)};
    wormIndex{nextDir(DIRECT)}     = 1:N(nextDir(DIRECT));

    tempcopy = array(wormIndex{:});

    if (obj.component == DIRECT); isVecNorm = -1; else; isVecNorm = 1; end;

    wormIndex{DIRECT} = amt:-1:1;
    wormIndex{nextDir(DIRECT)} = wormIndex{nextDir(DIRECT)}(end:-1:1);
end
    %--- Apply a constant shift ---%
    result = constant_shift(array,DIRECT,n,obj);
    
    %--- Fade the + edge of the array ---%
%    coeffSize               = ones(1,NSize);    coeffSize(DIRECT) = fadeSize;
%    coeff                   = zeros(coeffSize);
%    coeffIndex              = num2cell(ones(1,NSize));    coeffIndex{DIRECT} = 1:fadeSize;
%    coeff(coeffIndex{:})    = pchip([0 ceil(fadeSize/4) (fadeSize-1) fadeSize], [1 1 0 0], 1:fadeSize);
%    repMat                  = N;    repMat(DIRECT) = 1;
%    coeff                   = repmat(coeff,repMat);
%   
%    fadeIndex            = index;
%    fadeIndex{DIRECT}    = N(DIRECT);
%    fadeArray            = result(fadeIndex{:});
%    fadeArray            = repmat(fadeArray, coeffSize);
%    fadeIndex{DIRECT}    = (N(DIRECT)-fadeSize+1):N(DIRECT);
%    coeff                = flipdim(coeff,DIRECT);
%    result(fadeIndex{:}) = (1-coeff) .* result(fadeIndex{:}) ...
%                            + coeff .* fadeArray;

if dir == -1;
    result(wormIndex{:}) = isVecNorm*tempcopy;
end

end
