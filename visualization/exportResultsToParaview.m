function exportResultsToParaview(results, filename)
% This function will export array data (2D or 3D) from matlab into a format that ParaView can read.
% Specifically the file is of the type vts, which is a variant of the open source vtk format for
% data stored as a structured grid. The exported file is saved using the latest vtk xml binary form
% of the file specification format.
%
%>> <results>	imogen results data to be exported. Will find matching structure		struct
%					in base workspace if none is provided. 
%>> <filename>	file name to be written. Uses the name of the results argument if		str
%					none is provided.

	%-----------------------------------------------------------------------------------------------
	% Handle missing or improper input arguments
	%-------------------------------------------
	if (nargin < 1 || isempty(results) ) % Find results data is none provided
		vars	= evalin('base','who;');
		varTest = strfind(vars,'sx_');
		for i=1:length(varTest)
			if varTest{i};	results = evalin('base',sprintf('%s;',vars{i})); break; end
		end
		fprintf('WARNING: No data specified. Using %s.\n',vars{i});
	end
	
	if (nargin < 2 || isempty(filename) ); % Use argument name if none provided
		try
			filename = inputname(1); 
		catch ME
			filename = vars{i};
		end
		fprintf('WARNING: Using default file name of %s.\n',filename);
	end
	

	%-----------------------------------------------------------------------------------------------
	% Create the arrays needed for export
	%------------------------------------
	N = size(results.mass);
	if (length(N) < 3);		N(3) = 1; end
    ncells = numel(results.mass);
% 	offset = 0;
% 	comp = cell(1,5);
% 	compIndex = 1;
	comp = [];
	
	%-----------------------------------------------------------------------------------------------
	% Write to binary VTK XML file
	%-----------------------------
	fprintf('Exporting array data to VTK format:\n');
	fid = fopen([filename '.vtu'],'w','l','US-ASCII');
	generateVTKHeader(fid, N);
	
	fwrite(fid,sprintf('\n\t\t\t<CellData Scalars="mass" Vectors="velocity">'));
	
	
	%-----------------------------------------------------------------------------------------------
	% Write scalar arrays
	%--------------------
	fields = {'mass', 'mass';
			  'ener', 'energy';
			  'grav', 'gravity'};
	
	for i=1:length(fields)
		fprintf('\t%s..........',fields{i});
		
		if ~isfield(results,fields{i,1}) || isempty(results.(fields{i,1}))
			fprintf('empty.\n'); 
			continue
		end
		
		sca	 = reshape(results.(fields{i,1}),[1, ncells]);

% 		fwrite(fid,sprintf('\n\t\t\t\t<DataArray type="Float32" Name="%s" format="appended" offset="%g"/>', fields{i,2},offset));
% 		offset = offset + 4*length(sca);
% 		comp{compIndex} = sca; compIndex = compIndex + 1;
		
%		aStr = mat2str(sca,4);
		aStr = FastExportMat2Str(sca);
		aStr = aStr(2:(end-1));
		fwrite(fid,sprintf('\n\t\t\t\t<DataArray type="Float32" Name="%s" format="ascii">%s</DataArray>', fields{i,2},aStr));

        clear sca aStr;
		fprintf('complete.\n');
	end
	
	%-----------------------------------------------------------------------------------------------
	% Write vector arrays
	%--------------------
	fields = {'mom', 'velocity'; ...
			  'mag', 'magnetic_field'};
	comps = ['X', 'Y', 'Z'];
	
	for i=1:length(fields)
		fprintf('\t%s..........',fields{i});
		
		if ~isfield(results,[fields{i,1} comps(1)]) || isempty(results.([fields{i,1} comps(1)]))
			fprintf('empty.\n'); 
			continue
		end
		
		if strcmp(fields{i,1},'mom')
			for j=1:length(comps)
				results.([fields{i,1} comps(j)]) = results.([fields{i,1}  comps(j)]) ./ results.mass;
			end
		end
	
		vec = zeros(1,3*ncells);
		for j=1:3
			vec(j:3:end) = reshape(results.([fields{i,1} comps(j)]),[1, ncells]);
		end
		
% 		fwrite(fid,sprintf('\n\t\t\t\t<DataArray type="Float32" Name="%s" NumberOfComponents="3" format="appended" offset="%g"/>', fields{i,2},offset));
% 		offset = offset + 4*length(vec);
% 		comp{compIndex} = vec; compIndex = compIndex + 1;
		
%		vStr  = mat2str(vec,4);
		vStr  = FastExportMat2Str(vec);
		vStr  = vStr(2:(end-1));
		fwrite(fid,sprintf('\n\t\t\t\t<DataArray type="Float32" Name="%s" NumberOfComponents="3" format="ascii">%s</DataArray>', fields{i,2},vStr));

        clear vec vStr;
		fprintf('complete.\n');
	end
	
	fwrite(fid,sprintf('\n\t\t\t</CellData>'));
	
	
	fprintf('\n\tPositional data..........');
	grid = convertDGridToGrid(results.dGrid, size(results.mass));
	generateVTKPoints(fid, grid);
	fprintf('complete.\n');
	
	fprintf('\tCell data..........');
	generateVTKCells(fid, N);
	fprintf('complete.\n');
	
	generateVTKFooter(fid, comp);
	
	status = fclose(fid);
	switch status
		case -1;	fprintf('Unable to finalize file.\n');
		case 0;		fprintf('Export to vtk complete\n');
	end
		
end
