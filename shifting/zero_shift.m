function result = zero_shift(array, DIRECT, n, obj)
% Shifts an array so that the new cells added are identical to the pre-shifted edge cells.
%
%>> array		array to be shifted											double(Nx,Ny,Nz)
%>> DIRECT		direction/dimension along which to shift					int
%>> n			number of cells to shift									int
%>< obj			object owning the array										ImogenArray...
%<< result		shifted array												double(Nx,Ny,Nz)


    if isa(obj,'FluxArray')
    
        %--- Initialization ---%
        N		= size(array);
        NDim	= N(DIRECT);
        NSize	= length(N);
        index	= cell(1, NSize);
        for i=1:NSize,  index{i} = 1:N(i); end

        dir     = -sign(n); %Flip sign to correct inverse nature of shifting routines
        
        %--- SHIFT Array ---%
        if (dir>0)  % i = i + 1 (0->1, (N-1)->N)
            if strcmp(obj.fluxType, FluxArray.FLUXL)
                result = constant_shift(array, DIRECT, n, obj); return;
            else
                n        = abs(n);
                result   = zeros(N);
                inIndex  = index;   inIndex{DIRECT} = 1:(NDim-n);
                outIndex = index;   outIndex{DIRECT} = (n+1):NDim;
            end
        else        % i = i - 1 (1->0, (N+1)->N)
            if strcmp(obj.fluxType, FluxArray.FLUXR)
                result = constant_shift(array, DIRECT, n, obj); return;
            else
                n                         = abs(n);
                result                    = zeros(N);
                index{DIRECT}(1:(NDim-n)) = index{DIRECT}((n+1):NDim);
                inIndex  = index;   inIndex{DIRECT}  = (n+1):NDim;
                outIndex = index;   outIndex{DIRECT} = 1:(NDim-n);
            end
        end
        
        result(outIndex{:}) = array(inIndex{:}); return;
        
    else 
        result = constant_shift(array, DIRECT, n, obj); return;
    end
        
	
end