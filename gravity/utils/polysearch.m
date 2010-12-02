function X = polysearch(orders)

X = cell([numel(orders) 1]);

for q = 1:numel(orders)
    x0   = polyConditionerCoeff([8  8  8 ], 50, ones([1 orders(q)]));
    x0   = polyConditionerCoeff([12 12 12], 25 , x0);
    x0   = polyConditionerCoeff([16 16 16], 15 , x0);
    X{q} = polyConditionerCoeff([20 20 20], 5 , x0);
end

end
