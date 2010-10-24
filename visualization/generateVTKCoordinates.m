function compressedCoordinates = generateVTKCoordinates(fid, grid, compressed, offset)
% Generates the structured points portion of a vts file.
%
%>> arraySize			size of the array to be saved								double(3)
%>> dgrid				grid point locations for each dimension						cell(3)
%>> compressed			specifies if point data should be compressed				logical
%>> offset				# of bytes offset into the append array if compressed		int
%<< compressedPoints	compressed point data if appropriate						uint8

	if ( nargin < 4 || isempty(offset) );		offset = 0;				end
	if ( nargin < 3 || isempty(compressed) );	compressed = false;		end

	%--- Coordinate arrays (3D position values) ---%
	c = cell(1,3);
	[c{1}, c{2}, c{3}] = ndgrid(grid{1}, grid{2}, grid{3});
	npts = numel(c{1});
	
	fwrite(fid, sprintf('\t\t\t<Coordinates>'));
	if compressed;	compressedCoordinates = zeros(1,npts); end
	
	for i=1:3
		
		if compressed
			cStr = '';
			format = sprintf('format="appended" offset="%g"',offset);
			compressedPoints( ((i-1)*npts + 1):(i*npts) ) = dzip(points);
		else
			cStr = mat2str(reshape(c{i},[1, npts]),4);
			cStr = cStr(2:(end-1));
			format = 'format="ascii"';
			compressedPoints = [];
		end	
		
		fwrite(fid, sprintf('\n\t\t\t\t<DataArray type="Float32" %s>%s</DataArray>',format, cStr));

	end
	
	fwrite(fid, sprintf('\n\t\t\t</Coordinates>\n'));

	

end