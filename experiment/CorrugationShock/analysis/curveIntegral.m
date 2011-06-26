function result = curveIntegral(f,xi, dx)

result = 0;

for n = 1:2:(numel(f)-1)
    if (n+2) > numel(f)
        result = result + 3*(f(n+1)+f(n));
    else
        result = result + f(n) + 4*f(n+1) + f(n+2);
    end
end

result = result * dx / 3;

end