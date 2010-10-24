% A blurring diffusion-like operator
function rmat = simpleBlur(imat, coeff, n, specDir)
% This takes an input matrix imat and blurs it
%>> imat    The matrix to blur
%>> coeff   maps imat -> (1-coeff)*imat + coeff * blur
%>> n       Number of times to repeat the blur process
%>> specDir Specify direction - allows to blur in only 1 direction
%
%<< rmat    Result of the blurring

d = size(imat);
xset = [-1 0 1];
yset = [-1 0 1];
zset = [-1 0 1];

if nargin == 4;
    switch specDir
        case 1; yset = 0; zset = 0;
        case 2; xset = 0; zset = 0;
        case 3; xset = 0; yset = 0;
    end
end

if numel(d) < 3
    zset = [ 0 ];
    if d(2) == 1
         yset = [ 0 ];
    end
end

rmat = imat;
smat = zeros(size(imat));
q = 0;

for blurcount = 1:n
    for x = xset; for y = yset; for z = zset;
        if (x == 0) && (y == 0) && (z == 0) continue; end

        smat = smat + circshift(rmat, [x y z]) / sqrt(x^2 + y^2 + z^2);
        q = q + 1/sqrt(x^2 + y^2 + z^2);
    end; end; end;

    smat = smat / q; q = 0;

    rmat = rmat *(1-coeff) + smat * coeff;
    smat = zeros(size(imat));
end

end

