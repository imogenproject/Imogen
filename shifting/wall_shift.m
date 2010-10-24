function result = wall_shift(array, DIRECT, n, obj)
% Shifts an array so that the new cells added are constants defined by the initial conditions for 
% the array.
%
%>> array		array to be shifted											double((3),Nx,Ny,Nz)
%>> DIRECT		direction/dimension along which to shift					int
%>> n			number of cells to shift									int
%>< obj			object owning the array										ImogenArray...
%<< result		shifted array												double((3),Nx,Ny,Nz)

	%--- Initialization ---%
    result  = constant_shift(array, DIRECT, n, obj);
	N		= size(array);
    NDim	= N(DIRECT);
    NSize	= length(N);
    index	= cell(1, NSize);
    for i=1:NSize,  index{i} = 1:N(i); end
    
    dir                 = -sign(n); %Flip sign to correct inverse nature of shifting routines
    n                   = abs(n);
    repetitions         = ones(1,length(N));
    repetitions(DIRECT) = n;
    
	%--- Fill in wall ---%
	if (dir>0)
		index{DIRECT} = index{DIRECT}(1:n);
        edge = repmat(obj.transparentEdge(DIRECT,false), repetitions);
    else
		index{DIRECT} = index{DIRECT}((NDim-n+1):NDim);
        edge = repmat(obj.transparentEdge(DIRECT,true), repetitions);
	end
	
	result(index{:}) = edge;
end
