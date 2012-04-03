
x = 0:.01:1;
y = exp(-3*x);

[kre kim] = ndgrid(0:.04:2.5, 0:.04:5);

f = zeros(size(kre));
g = zeros(size(f));
h = zeros(size(g));


parfor n = 1:numel(f);
    [f(n) g(n) h(n)] = multiwaveFit(y, .01, kre(n) + 1i*kim(n), 1);
end


