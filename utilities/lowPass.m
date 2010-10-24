function lp = lowPass(f, c)

lp(1) = f(1);

for t = 2:numel(f)
    lp(t) = lp(t-1) + c*(f(t) - lp(t-1));
end

lp = reshape(lp, size(f));

end