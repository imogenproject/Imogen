function result = mirror_shift(array, DIRECT, n, obj)
% Shifts an array reflectively, so that the cells lost at one end are added as the new cells gained 
% after being reflected across the shifting axis.
%
%>> array		array to be shifted	                                    	double((3),Nx,Ny,Nz)
%>> DIRECT		direction/dimension along which to shift	            	int
%>> n	    	number of cells to shift	                            	int
%>< obj	    	object owning the array	                                	ImogenArray...
%<< result		shifted array	                                        	double((3),Nx,Ny,Nz)


    %--- Initialization ---%
	N	    = size(array);
    NDim	= N(DIRECT);
    NSize	= length(N);
    index	= cell(1, NSize);
    for i=1:NSize,  index{i} = 1:N(i); end
    
    dir = -sign(n); %Flip sign to correct inverse nature of shifting routines
    n	= abs(n);

    %--- SHIFT Array ---%
    if (dir>0)
    	index{DIRECT}((n+1):NDim)       = index{DIRECT}(1:(NDim-n));
        flipIndex                       = 1:n;
    	index{DIRECT}(flipIndex)        = (n+1):-1:2;
	else
    	index{DIRECT}(1:(NDim-n))     	= index{DIRECT}((n+1):NDim);
        flipIndex                       = (NDim-n+1):NDim;
    	index{DIRECT}(flipIndex)        = (NDim-1):-1:(NDim-n);
    end
    
    result = array(index{:});
    
    %--- Flip sign on appropriate vector indeces ---%
    if (obj.component == DIRECT)
        index{DIRECT} = flipIndex;
        result(index{:}) = -result(index{:});
    end

end