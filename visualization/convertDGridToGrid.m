function grid = convertDGridToGrid(dgrid, arraySize)
% This function converts an Imogen grid spacing cell array into a grid cell array that contain 
% absolute positions for each spatal direction.
%
%>> dgrid		grid spacing array													cell(3)
%>> arraySize	# of grid points (not cells) in each spatial dimension				double(3)
%<< grid		grid point positions												cell(3)

	if length(arraySize) < 3;	arraySize = [arraySize 1];	end
	if isa(dgrid,'double');		dgrid = num2cell(dgrid);	end
	
	grid = cell(1,3);
	maxVal = 0;
	for i=1:3
		grid{i} = zeros(1,arraySize(i)+1);
		if length(dgrid{i}) == 1
			grid{i}(2:end)  = dgrid{i};
			grid{i}			= cumsum(grid{i});
		else
			fills = num2cell(ones(1,3)); fills{i} = 1:arraySize(i);
			grid{i}(2:end) = cumsum(dgrid{i}(fills{:}));
		end
		maxVal = max(maxVal, max(grid{i}));
	end
	
	% Rescale so that largest dimension is between 0 and 100
	for i=1:3
		grid{i} = 100/maxVal .* grid{i};
	end
	
	

end