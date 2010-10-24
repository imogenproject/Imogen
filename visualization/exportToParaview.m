function exportToParaview(arrays, dGrid, filename, tags)
% This function will export array data (2D or 3D) from matlab into a format that ParaView can read.
% Specifically the file is of the type vts, which is a variant of the open source vtk format for
% data stored as a structured grid. The exported file is saved using the latest vtk xml binary form
% of the file specification format.
%
%>> arrays		data arrays to be exported. Can be either a single array or a cell		double/cell
%				array of arrays to export.
%>> <dGrid>		grid spacings between cells. Defaults to [1 1 1].						double
%>> <filename>	file name to be written. A name is generated if none is provided.		str
%>> <tag>		tags to use for the array(s) in the vtk file. If none is provided		str/cell
%				one will be generated.

	%-----------------------------------------------------------------------------------------------
	% Handle missing or improper input arguments
	%-------------------------------------------	
	if ~iscell(arrays),		arrays = {arrays}; end
    if ~iscell(tags),       tags   = {tags};   end
	
	if (nargin < 4 || isempty(tags));
		tags = cell(1,length(arrays));
		for i=1:length(arrays)
			tags{i} = inputname(1) + '_' + i;
		end
		fprintf('WARNING: Tags have been auto-generated.\n');
	end

	if (nargin < 3 || isempty(filename) ); % Use argument name if none provided
		filename = inputname(1); 
		fprintf('WARNING: Using default file name of %s.\n',filename);
	end
	
	if (nargin < 2 || isempty(dGrid) );
		dGrid = [1 1 1];
		fprintf('WARNING: Using default grid spacing of [1 1 1].\n');
	end
	
	%-----------------------------------------------------------------------------------------------
	% Write to binary VTK XML file
	%-----------------------------
	fprintf('Exporting array data to VTK format:\n');
	fid = fopen([filename '.vtu'],'w','l','US-ASCII');
	
	N = size(arrays{1});
	if (length(N) > 3 || N(1) == 3), N = N(2:end); end %Remove vector components from array size.
	generateVTKHeader(fid, N);
	fwrite(fid,sprintf('\n\t\t\t<CellData>'));
	
	%-----------------------------------------------------------------------------------------------
	% Write array
	%------------
	for i=1:length(arrays)
		fprintf('\t%s..........',tags{i});
		
		N			= size(arrays{1});
		ncells		= numel(arrays{1});
		isaVector	= (N(1) == 3);
		
		%--- Make 2D arrays into 3D arrays ---%
		if (~isaVector && length(N) < 3),		N(3) = 1; end
		if (isaVector && length(N) < 4),		N(4) = 1; end
		
		%-------------------------------------------------------------------------------------------
		% Write scalar arrays
		%--------------------
		if ~isaVector

			sca	 = reshape(arrays{i},[1, ncells]);
			aStr = mat2str(sca,4);
			aStr = aStr(2:(end-1));
			fwrite(fid,sprintf('\n\t\t\t\t<DataArray type="Float32" Name="%s" format="ascii">%s</DataArray>', tags{i}, aStr));
			clear sca aStr;
			fprintf('complete.\n');

		%-------------------------------------------------------------------------------------------
		% Write vector arrays
		%--------------------
		else
			vec = zeros(1,3*ncells);
			for j=1:3
				vec(j:3:end) = reshape(squeeze(arrays{i}(j,:)),[1, ncells]);
			end
			vStr  = mat2str(vec,4);
			vStr  = vStr(2:(end-1));
			fwrite(fid,sprintf('\n\t\t\t\t<DataArray type="Float32" Name="%s" NumberOfComponents="3" format="ascii">%s</DataArray>', tags{i}, vStr));

			clear vec vStr;
			fprintf('complete.\n');
		end
	end
	
	
	fwrite(fid,sprintf('\n\t\t\t</CellData>'));
	
	fprintf('\n\tPositional data..........');
	grid = convertDGridToGrid(dGrid, size(arrays{1}));
	generateVTKPoints(fid, grid);
	fprintf('complete.\n');
	
	fprintf('\tCell data..........');
	generateVTKCells(fid, N);
	fprintf('complete.\n');
	
	comp = [];
	generateVTKFooter(fid, comp);
	
	status = fclose(fid);
	switch status
		case -1;	fprintf('Unable to finalize file.\n');
		case 0;		fprintf('Export to vtk complete\n');
	end
		
end
