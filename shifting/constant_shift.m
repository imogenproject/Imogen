function result = constant_shift(array, DIRECT, n, obj)
% Shifts an array so that the new cells added are identical to the pre-shifted edge cells.
%
%>> array        array to be shifted                                 double((3),Nx,Ny,Nz)
%>> DIRECT       direction/dimension along which to shift            int
%>> n            number of cells to shift                            int
%>< obj          object owning the array                             ImogenArray...
%<< result       shifted array                                       double((3),Nx,Ny,Nz)

    %--- Initialization ---%
    N            = size(array);
    NDim         = N(DIRECT);
    NSize        = length(N);
    index        = cell(1, NSize);
    for i=1:NSize,  index{i} = 1:N(i); end
    
    dir = -sign(n); %Flip sign to correct inverse nature of shifting routines
    n        = abs(n);

        %--- SHIFT Array ---%
        if (dir>0)
                index{DIRECT}((n+1):NDim) = index{DIRECT}(1:(NDim-n));
                index{DIRECT}(1:n)        = 1;
        else
                index{DIRECT}(1:(NDim-n))      = index{DIRECT}((n+1):NDim);
                index{DIRECT}((NDim-n+1):NDim) = NDim;
        end
        
        result = array(index{:}); return;
end
