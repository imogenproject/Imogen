function exportCompressedResultsToParaview(results, filename)
% This function will export array data (2D or 3D) from matlab into a format that ParaView can read.
% Specifically the file is of the type vts, which is a variant of the open source vtk format for
% data stored as a structured grid. The exported file is saved using the latest vtk xml binary form
% of the file specification format.
%
%>> results		results data to be exported											struct
%>> filename	name of the file to written											str

	%-----------------------------------------------------------------------------------------------
	% Create the arrays needed for export
	%------------------------------------
	N = size(results.mass);
	if (length(N) < 3),	N(3) = 1; end
    npts   = prod(N);
	offset = 0;
	compData = cell(1,10);
	compIndex = 1;
	
	%-----------------------------------------------------------------------------------------------
	% Write to binary VTK XML file
	%-----------------------------
	fprintf('Exporting array data to VTK format.......');
	fid = fopen([filename '.vts'],'w');
	generateVTKHeader(fid, N);
	
	[pStr, compData{compIndex}] =  generateVTKPoints(N,true,offset);
	compIndex = compIndex + 1;
	fwrite(fid,pStr);
	offset = offset + 3*npts*4;
	
	fwrite(fid,sprintf('\t\t\t<PointData Scalars="mass" Vectors="velocity">\n'));
	
	%-----------------------------------------------------------------------------------------------
	% Write scalar arrays
	%--------------------
	fields = {'mass', 'mass'; ...
			  'ener', 'energy'};
	
	for i=1:length(fields)

		compData{compIndex}	= dzip(reshape(results.(fields{i,1}),[1, npts]));
		compIndex = compIndex + 1;

		fwrite(fid,sprintf(['\t\t\t\t<DataArray type="Float32" Name="%s" format="appended" ' ...
						   'offset="%g"></DataArray>\n'], fields{i,2},offset));
		offset = offset + npts*4;
	end
	
	%-----------------------------------------------------------------------------------------------
	% Write vector arrays
	%--------------------
	fields = {'mom', 'velocity'; ...
			  'mag', 'magnetic_field'};
	comps = ['X', 'Y', 'Z'];
	
	for i=1:length(fields)
		
		if strcmp(fields{i,1},'mom')
			for j=1:length(comps)
				results.([fields{i,1} comps(j)]) = results.([fields{i,1}  comps(j)]) ./ results.mass;
			end
		elseif strcmp(fields{i,1},'mag')
			if isempty(results.([fields{i,1} comps(1)])), continue; end
		end
	
		vec = zeros(1,3*npts);
		for j=1:3
			vec(j:3:end) = reshape(results.([fields{i,1} comps(j)]),[1, npts]);
		end

		compData{compIndex} = dzip(vec);
		compIndex = compIndex + 1;
		
		fwrite(fid,sprintf(['\t\t\t\t<DataArray type="Float32" Name="%s" NumberOfComponents="3" ' ...
			'format="appended" offset="%g"></DataArray>\n'], fields{i,2},offset));
		
		offset = offset + 3*npts*4;
		
	end
	
	fwrite(fid,sprintf('\t\t\t</PointData>\n'));
	generateVTKFooter(fid,compData);
	
	status = fclose(fid);
	switch status
		case -1;	fprintf('Unable to finalize file.\n');
		case 0;		fprintf('Export to vtk complete\n');
	end
		
end