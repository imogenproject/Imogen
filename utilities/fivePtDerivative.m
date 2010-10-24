function result = fivePtDerivative(imoArray,dir,dGrid)
%   Computes the five point centered derivative of the array in the specified direction. The five 
%   point derivative is: (A(x-2dx) - 8A(x-dx) + 8A(x+dx) - A(x+2dx))/(12dx)
%
%>< imoArray	array to compute the derivative on							ImogenArray			H
%>> dir			index direction of the derivative							int
%>> dGrid		grid separation vector										cell(3)
%<< result		derivative of the imoArray.array							double(?)

    %--- Compute the derivative ---%
    result = ( -imoArray.shift(dir, 2)    + 8*imoArray.shift(dir, 1) ...
               -8*imoArray.shift(dir, -1) + imoArray.shift(dir, -2) ) / (12*dGrid{dir});
end