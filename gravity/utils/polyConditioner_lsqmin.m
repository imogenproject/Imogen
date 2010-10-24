function lambda = polyConditioner_lsqmin(order, iterMax)
% minimize del_{deltaset}

deltaset = ones([1 order]);
h = .001;

partials = zeros([1 order]);
niters = 0;

lamint = @(deltas) quad(@(x) (1 - q(x, deltas)).^2, 0, 2);
Z = [];


while niters < iterMax
    lambdaIntegral0 = lamint(deltaset);
    for z = 1:order
        deltaset(z) = deltaset(z) + h;
        partials(z) = (lamint(deltaset) - lambdaIntegral0)/sqrt(h);
        deltaset(z) = deltaset(z) - h;
    end
    
    for z = 1:order
        deltaset(z) = deltaset(z) - sign(partials(z))*min(abs([3*h sqrt(abs(partials(z)))]) ); 
    end
niters = niters + 1;
Z(niters,:) = deltaset(:);

end

plot(Z);
disp(deltaset);    

end

function F = q(lambda, deltaset)

F=zeros(size(lambda));

for t = 1:numel(deltaset)
   F = F + deltaset(t)*lambda.^t;
end

end