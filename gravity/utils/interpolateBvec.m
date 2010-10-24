function result = interpolatedBvec(rhos, poss, bvec0, cconst)

q = numel(poss);
h = poss{q}(1,2,1,1) - poss{q}(1,1,1,1);

dims = bvec0(7:9)+1;

bvec = bvec0;
trimme = [0 0 0];

isfirstnonsingle = 1;
subme = [0 0];

for q = 1:3;
    if mod(dims(q),2) == 0 % even
        bvec(6+q) = dims(q)/2;
        bvec(3+q) = bvec(3+q) + h;
        trimme(q) = 1;
        if dims(q) > 1;
            if isfirstnonsingle; subme(1) = 1; isfirstnonsingle = 0; else subme(2) = 1; end;
        end
    else % odd
        bvec(6+q) = (dims(q)-1)/2;
    end

    if bvec0(6+q) == 0; bvec(6+q) = 0; trimme(q) = 0; end % We only want to make a single plane

    if bvec0(6+q) ~= 0
        bvec(q) = bvec(q)-2*h;
        bvec(3+q) = bvec(3+q) + 2*h;
        bvec(6+q) = bvec(6+q) + 2;
    end
end

% Generate values to interpolate using.
phi0 = mg_bc_matlab(rhos, poss, bvec, cconst);
phi0 = squeeze(phi0);

phipr = interp2(phi0, 'cubic');

zeta = size(phipr);
zeta = zeta - 2 - subme;

result = reshape(phipr(3:zeta(1), 3:zeta(2)), dims);

end
