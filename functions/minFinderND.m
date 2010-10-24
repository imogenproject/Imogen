function [minVal, minVecIndex] = minFinderND(array)
%   Finds the minimum value of an N dimensional arrray and the corresponding index of the first
%   array dimension. The index for the minimum value along the first array dimension can be used
%   to determine X, Y or Z directions for vector arrays.
%
%>> array		array for the minimum value search									double
%<< minVal        scalar min value for the entire array								double
%<| minVecIndex   index of the first array dimension where min occurs				int


    %--- Find array dimension ---%
    DIM = max(size(size(array)));
    
    %--- Search for maximum by dimensions ---%
    for i=0:(DIM-2), array = min(array,[],DIM-i); end % Iterate from DIM to 1 so no need for squeeing.
    [minVal, minVecIndex] = min(array,[],1); % Find the min of the last remaining (first) dimension

end