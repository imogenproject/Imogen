function generateVTKPoints(fid, grid)
% Generates the unstructured points portion of a vtu file.
%
%>> fid			file identifier for writing											handle
%>> grid			grid coordinate positions										cell(3)

	%--- Points array (3D position values) ---%
	p = cell(1,3);
    [p{1}, p{2}, p{3}] = ndgrid(grid{1}, grid{2}, grid{3});
	npts = numel(p{1});
	
    points = zeros(1,3*npts);
    points(1:3:end) = reshape(p{1},[1, npts]);
    points(2:3:end) = reshape(p{2},[1, npts]);
    points(3:3:end) = reshape(p{3},[1, npts]);
    clear('p');
	
%	pStr = mat2str([points(1:3) points],4);
	pStr = FastExportMat2Str([points(1:3) points]);
	pStr = pStr(2:(end-1));
	
	fwrite(fid, sprintf(['\n\t\t\t<Points>\n\t\t\t\t<DataArray type="Float32" NumberOfComponents="3"', ... 
					 blanks(1), 'format="ascii">%s</DataArray>\n\t\t\t</Points>'],pStr) ...
			 );
end
