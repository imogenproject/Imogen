function result = transparent_shift(array, DIRECT, n, obj)
% Shifts an array so that the new cells are spline interpolated values designed for transparency.
%
%>> array       array to be shifted.											double(GRID)
%>> DIRECT		direction/dimension along which to shift                        int
%>> n			number of cells to shift                                        int
%>< obj			object owning the array                                         ImogenArray     H
%<< result		shifted array                                                   double(GRID)

	%--- Initialization ---%
	N		= size(array);
    NDim	= N(DIRECT);
    NSize	= length(N);
    index	= cell(1, NSize);
    for i=1:NSize,  index{i} = 1:N(i); end
    
    dir = -sign(n); %Flip sign to correct inverse nature of shifting routines
    n	= abs(n);
	permShift = { [2 3 1], [3 1 2], [1 2 3] };
	
	
	%--- SHIFT Array ---%
	if (dir>0)	
		
		% Shift bulk of array by n such that i' = i-1
		index{DIRECT}((n+1):NDim)	= index{DIRECT}(1:(NDim-n));
		index{DIRECT}(1:n)			= 1;
		result = array(index{:});

		% Extract the edge from edge from the array for interpolation
		inSize			= N;
		inSize(DIRECT)	= 5;
		inVals = zeros(inSize); 
		inIndex			= index; 
		inIndex{DIRECT} = 1:3;																		
		index{DIRECT}	= min( (3:-1:1)+n, NDim );
		if ~isa(array,'double'),			inVals(inIndex{:}) = gather( result(index{:}) ); %r2009b: iscodistributed
		else								inVals(inIndex{:}) = result(index{:});
		end

		% If storing edge values use them instead of zero.
		if (obj.edgeStore.ACTIVE(DIRECT))
			inIndex{DIRECT}		 = 4; 
			index{DIRECT}		 = 5;
			inVals(inIndex{:})	 = obj.transparentEdge(DIRECT,false);
			inVals(index{:})	 = inVals(inIndex{:});
		end

		% Interpolate and populate result array
		inVals	= permute(inVals, permShift{DIRECT});
		outVals = pchip( [1:3 obj.bcInfinity (obj.bcInfinity+1)], inVals, 3+(n:-1:1) );
		outVals = ipermute(outVals, permShift{DIRECT});
		
		index{DIRECT} = 1:n;
		
	else
		% Shift bulk of array by n such that i' = i+1
		index{DIRECT}(1:(NDim-n))		= index{DIRECT}((n+1):NDim);
		index{DIRECT}((NDim-n+1):NDim)	= NDim;
		result = array(index{:});

		% Extract the edge from edge from the array for interpolation
		inSize			= N; 
		inSize(DIRECT)	= 5;
		inVals = zeros(inSize); 
		inIndex			= index; 
		inIndex{DIRECT} = 1:3; 
		index{DIRECT}	= max( NDim-(2:-1:0)-n, 1 );
		if ~isa(array,'double'),			inVals(inIndex{:}) = gather( result(index{:}) ); %r2009b: iscodistributed
		else								inVals(inIndex{:}) = result(index{:});
		end

		% If storing edge values use them instead of zero.
		if (obj.edgeStore.ACTIVE(DIRECT))
			inIndex{DIRECT}		= 4; 
			index{DIRECT}		= 5;
			inVals(inIndex{:})	= obj.transparentEdge(DIRECT,true);
			inVals(index{:})	= inVals(inIndex{:});
		end

		% Interpolate and populate result array
		inVals  = permute(inVals, permShift{DIRECT});
		outVals = pchip( [1:3 obj.bcInfinity (obj.bcInfinity+1)], inVals, 3+(1:n) );          
		outVals = ipermute(outVals, permShift{DIRECT});

		index{DIRECT} = (NDim-n+1):NDim;
		
	end
	
	%--- Ignore small value cases ---%
    try
        test = abs(array(index{:}) - outVals) < 1e-4;
    catch MERR
        outVals = outVals';
        test = abs(array(index{:}) - outVals) < 1e-4;
    end
	outVals = not(test).*outVals + test.*array(index{:});
	
	result(index{:}) = outVals;
	
end
