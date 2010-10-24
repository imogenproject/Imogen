function X = polysearch(orders)

X = cell([numel(orders) 1]);

for q = 1:numel(orders)
    x0   = polyConditionerCoeff([8  8  8 ], 100, ones([1 orders(q)]));
    x0   = polyConditionerCoeff([12 12 12], 50 , x0);
    X{q} = polyConditionerCoeff([16 16 16], 25 , x0);
end

end
