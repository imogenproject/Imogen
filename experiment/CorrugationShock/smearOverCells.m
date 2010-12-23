function result = smearOverCells(data, n)

x0 = data(1,1,1);
x1 = data(end,1,1);

q = 1:size(data,1);
qh = numel(q)/2;

n = n/5; % this looks good

y = x0 + (x1-x0)*exp((q - qh)/n) ./ (1 + exp((q-qh)/n));

result = data;
for x = 1:numel(q)
	result(x,:,:) = y(x);
end

end
