function C = polyCondMatrix(x1, x2, x3, x4)
% This is a bit of a bastard... You can edit it easily to poop out the results from searching
% a gridded volume of parameter space for any number of parameters.

[a b c d] = ndgrid(x1, x2, x3, x4);

C = zeros(size(c));

parfor q = 1:numel(C);
    C(q) = polyCondCheck([16 16 16], [a(q) b(q) c(q) d(q)]);
end

end
