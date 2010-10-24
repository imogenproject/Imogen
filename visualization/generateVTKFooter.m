function generateVTKFooter(fid, compData)
% Generates and write the footer portion of a vtk file.
%
%>> fid			file identifier for writing											handle
%>> compData	data arrays to be written to the appended data section				cell(?)

	if (nargin < 2);	compData = []; end
	
	if isempty(compData)
		fwrite(fid,sprintf('\n\t\t</Piece>\n\t</UnstructuredGrid>\n</VTKFile>'));
	else
		
		fwrite(fid,sprintf('\n\t\t</Piece>\n\t</UnstructuredGrid>\n\t<AppendedData encoding="raw">_'));
		
		for i=1:length(compData)
			if isempty(compData{i}); continue; end
			fwrite(fid, length(compData{i})*4,'uint');
			fwrite(fid, compData{i},'float32');
		end
		
		fwrite(fid,sprintf('</AppendedData>\n</VTKFile>'));
		
	end
		
		

end