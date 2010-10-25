function sumVal = sumND(array)
% Finds the total sum of an N dimensional arrray.
%
%>> array     array to search through for a maximum value					double  [n1 n2 ...]
%<< sumVal    sum value for the entire array								double  #


    %--- Find array dimension ---%
    DIM = length(size(array));
    
    %--- Search for maximum by dimensions ---%
    for i=0:(DIM-2), array = sum(array,DIM-i); end % Iterate from DIM to 1 so no need for squeeing.
    sumVal = sum(array,1);
    
end