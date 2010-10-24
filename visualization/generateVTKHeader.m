function generateVTKHeader(fid, arraySize)
% Generates the header portion of an unstructured grid vtk file.
%
%>> fid			file identifier for writing											handle
%>> arraySize	size of the arrays to be written									int(2 or 3)
		
	Nc = prod(arraySize);
	if (length(arraySize) < 3); arraySize(3) = 1; end
	Np = prod(arraySize+1);

	fwrite(fid, sprintf('<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">') );
	fwrite(fid, sprintf('\n\t<UnstructuredGrid>\n\t\t<Piece	NumberOfPoints="%g" NumberOfCells="%g">', Np+1, Nc));
	
end