function [maxVal, maxVecIndex] = maxFinderND(array)
%   Finds the maximum value of an N dimensional arrray and the corresponding index of the first
%   array dimension. The index for the maximum value along the first array dimension can be used
%   to determine X, Y or Z directions for vector arrays.
%
%>> array			array for maximum value	search								double
%<< maxVal			maximum value for the searched array						double
%<| maxVecIndex		index of the first array dimension where maxVal occurs		int

    if isa(array,'GPUdouble') == 1
     [maxVal, maxVecIndex] = maxFinderND(double(array));
      return;
    end

    %--- Find array dimension ---%
    DIM = ndims(array);
    
    %--- Search for maximum by dimensions ---%    

	for i=0:(DIM-2), array = max(array,[],DIM-i); end

	if ~(isa(array,'double') ||  isa(array,'uint8')); array = gather(array); end %r2009b: iscodistributed

	[maxVal, maxVecIndex] = max(array,[],1);
end
