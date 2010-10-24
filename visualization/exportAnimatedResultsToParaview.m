function exportAnimatedResultsToParaview(baseFileName, padDigits, varargin)
% Sequentially exports the varargin list of results structures to the Paraview VTK XML format for 
% animated data analysis through sequential numbering of the export files. The provided file name
% is appended with numeric identifiers so that filename.vts becomes filename_###.vts uniquely for 
% each file exported.
%
%>> baseFileName	full path to the base file name to use during export				str
%>> padDigits       Number of digits to use in padding the export numbers.              int
%>> varargin		results structures to be exported									struct OR
%>>                 if varargin is a single double entry, it is taken to be the offset  double
	
    %--- Initialization ---%
    fprintf('Beginning file export...\n');
    vlen = length(varargin);

	%--- Find variables ---%
	%		If variables aren't specified in input argument then search and import them from the 
	%		base workspace.    
	if ( nargin < 3 || isempty(varargin) || ~ischar(varargin{1}) )
        fprintf('No files specified. Importing from workspace.\n');		
		buffer = evalin('base','who(''sx_*'')');	buffer = sort(buffer);
		vlen = length(buffer);
		
		index = 1;
		vars = cell(1,vlen);
		for i=1:vlen
			startRes = strfind(buffer{i},'START');
			if ~isempty(startRes)
				vars{1} = evalin('base',buffer{i}); 
				index = 2; 
			end
			
			finalRes = strfind(buffer{i},'FINAL');
			if ~isempty(finalRes)
				vars{end} = evalin('base',buffer{i}); 
			end
		end
		
		for i=1:vlen
			startRes = strfind(buffer{i},'START');
			finalRes = strfind(buffer{i},'FINAL');
			if isempty(startRes) && isempty(finalRes)
				vars{index} = evalin('base',buffer{i});
				index = index + 1;
			end
        end
        
        if ~isempty(varargin) && isa(varargin{1},'double');     offset = varargin{1};
        else                                                    offset = 0;
        end
        
	else
		vars = varargin;
        offset = 0;
    end

    %--- Determine index padding ---%
    if (nargin < 2 || isempty(padDigits))
        temp = vlen+offset;
        for i=1:10
            temp = floor(temp / 10);
            if (temp == 0), pLen = i + 1; break; end
        end
    else
        pLen = padDigits;
    end
	
  	%--- Determine save path ---%
	if (nargin < 1 || isempty(baseFileName))
        pathstr = '';
        name    = 'Imogen_Data';
	else
        [pathstr, name] = fileparts(baseFileName);
	end
    
	%--- Export sequential files ---%
	for i=1:vlen
        index = i+offset-1;
		filename = [pathstr, name, '_', Paths.paddedNumber(index,pLen)];
		varname = num2str(vars{i}.iter);
		fprintf('(%g/%g) Exporting %s to %s:  ',i,vlen,varname,filename);
		exportLegacyVTK(vars{i}, filename);
%		exportResultsToParaview(vars{i},filename);
	end

	fprintf('Export complete.\n');
	
end
