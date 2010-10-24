function generateVTKCells(fid, arraySize)
% Generates the cell data (connectivity, offsets, and cell types) for unstructured grids.
%
%>> fid			file identifier for writing											handle
%>> arraySize	size of the arrays to be written									int(2 or 3)


	%--- Initialization ---%
	Nc	    = arraySize;		ncells  = prod(Nc);
	Np		= Nc+1;				npts	= prod(Np);
	cellIndex = 1;
	c = zeros(1,8*ncells);
	
	%--- Generate connectivity array ---%
	for k=1:Nc(3)
		for j=1:Nc(2)
			for i=1:Nc(1)
				cmin = 8*(cellIndex-1) + 1;
				cmax = 8*cellIndex;
				xs = i*ones(1,8) + [0 1 1 0 0 1 1 0];
				ys = j*ones(1,8) + [0 0 1 1 0 0 1 1];
				zs = k*ones(1,8) + [0 0 0 0 1 1 1 1];
				c(cmin:cmax) = xs + (ys-1).*Np(1) + (zs-1).*Np(1).*Np(2);
				cellIndex = cellIndex + 1;
			end
		end
	end
			

	%--- Write cell data arrays --%
	fwrite(fid, sprintf('\n\t\t\t<Cells>'));
	
%	connects = mat2str( c );
	connects = FastExportMat2Str_int(c);
	connects = connects(2:(end-1));
	fwrite(fid, sprintf('\n\t\t\t\t<DataArray type="Int32" Name="connectivity" format="ascii">%s</DataArray>', connects) );
	clear c connects;
    
%	offsets = mat2str( 8*(1:ncells) );
	offsets = FastExportMat2Str_int( 8*(1:ncells) );
	offsets = offsets(2:(end-1));
	fwrite(fid, sprintf('\n\t\t\t\t<DataArray type="Int32" Name="offsets" format="ascii">%s</DataArray>', offsets) );
	clear offsets;
    
%	types = mat2str( 12*ones(1,ncells) );
	types = FastExportMat2Str_int( 12* ones(1,ncells) );
	types = types(2:(end-1));
	fwrite(fid, sprintf('\n\t\t\t\t<DataArray type="UInt8" Name="types" format="ascii">%s</DataArray>', types) );
	clear types;
    
	fwrite(fid, sprintf('\n\t\t\t</Cells>'));

end
