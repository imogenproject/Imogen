classdef ImogenRecord < handle
% Class for creating data records of data. Includes a number of routines for converting objects
% into strings for writing to file or viewing in the command panel.
	
%===================================================================================================
	properties (SetAccess = private, GetAccess = private, Transient = true) %	P R I V A T E    [P]

	end%PRIVATE
	
%===================================================================================================
    properties (SetAccess = public, GetAccess = public, Transient = true) %			P U B L I C  [P]
		source;	
	end%PUBLIC
	
%===================================================================================================
    methods (Access = public) %														P U B L I C  [M]
		
		function obj = ImogenRecord(src)
			if ( nargin == 1 && ~isempty(src) )
				obj.source = ImogenRecord.valueToString(src); 
			end
		end
		
		
		
	end%PUBLIC
	
%===================================================================================================
	methods %																	G E T / S E T	 [M]
	
		function set.source(obj,src)
			obj.source = ImogenRecord.valueToString(src);
		end
		
		function res = get.source(obj)
			res = obj.source;
		end
		
	end%GET/SET
	
%===================================================================================================
    methods (Access = private) %												P R I V A T E    [M]
		
	end%PRIVATE
	
%===================================================================================================
	methods (Static = true) %														S T A T I C  [M]
%___________________________________________________________________________________________________ valueToString
% Converts a non-string value to string.
%>>	inVal	non-string input value														*
%>>	<skips>	array of field names to skip												cell
%>> <depth> tab depth for the string lines when reported								int
%<<	value	string resulted from conversion of the inVal.								str
		function value = valueToString(inVal, skips, depth)
			if (nargin < 2 || isempty(skips)); skips = {''}; end
			
			switch (class(inVal))
				case 'struct';				value = ImogenRecord.structToStr(inVal);
					
				case 'char';				value = inVal;
					
				case 'cell';				value = ImogenRecord.cellToStr(inVal);
					
				case 'logical';				value = ImogenRecord.logicalToStr(inVal);
					
				case {'double','single','uint8'};
											value = ImogenRecord.numberToStr(inVal);
											
				case 'codistributed';		value = ImogenRecord.distToStr(inVal);
					
				case 'function_handle';		value = func2str(inVal);
					
				case 'codistributor1d';		value = ImogenRecord.codistToStr(inVal);
					
				otherwise %Assumed to be a class object
					if isobject(inVal)
						if (nargin < 3 || isempty(depth)); depth = 0; end
							objLen = length(inVal);
							value = '';
							if objLen == 1
								value = ImogenRecord.objectToStr(inVal, skips, depth);
							else
								for i=1:length(inVal)
									value = [value, sprintf('\n\n---++++ Index: %g\n',i)];
									value = [value, ImogenRecord.objectToStr(inVal(i),skips, depth)];
								end
							end
					else
						warning('ImogenRecord:UknownValueType','Unknown value type supplied.');
						value = '';
					end
			end
		end
		
%___________________________________________________________________________________________________ numberToStr
% Converts a numeric type into a string.
% * inNum	the input number or number array											[numeric]
% # value	the string resulting from conversion of the number or matrix array.			str
		function value = numberToStr(inNum)
			if ( (ndims(inNum) < 3) && (numel(inNum) < 13) )
				value = mat2str(inNum);
				info = '';
			else
				value = strcat('Matrix ',mat2str(size(inNum)) );
				info =	[' | min: ', num2str(minFinderND(inNum)), ...
						' | max: ', num2str(maxFinderND(inNum))];
			end
			value = ['(',class(inNum),') ', value,	info];
		end
		
%___________________________________________________________________________________________________ distToStr
% Converts a codistributed type into a string.
% * inDist	the input matrix or number array											[numeric]
% # value	the string resulting from conversion of the number or matrix array.			str
		function value = distToStr(inDist)
			value = ['size: ' mat2str(size(inDist))];
			value = ['(',class(inDist),') ', value];
		end	
		
%___________________________________________________________________________________________________ codistToStr
% Converts a codistributor type into a string
% * inCodist	the input codistributor1d object										codistributor1d
% # value	the string resulting from conversion of the codistributor.					str
		function value = codistToStr(inCodist)
			value = class(inCodist);
		end
		
%___________________________________________________________________________________________________ objectToStr
% Converts an object instance into a string.
% * inObj	the input object															Object
% # value	the string resulting from conversion of the cell array.						Str
		function value = objectToStr(inObj, skips, depth)
			depth = depth+1;
			prefix = [blanks(3*depth) '* '];
			names = fieldnames(inObj);
			empty = sprintf('%sEMPTY:',prefix);
			value = sprintf('\n%sCLASS: %s\n', prefix, class(inObj));
			for i=1:length(names)
				%Skip constants and specified fields
				if (lower(names{i}(1)) ~= names{i}(1)) || any(strcmpi(names{i}, skips)); continue; end
				
				val = inObj.(names{i});
				if isempty(val)
					empty = strcat(empty, ' | ', names{i}); continue;
				else
					value = [value, sprintf('%s%s: %s\n', prefix, names{i}, ...
							 ImogenRecord.valueToString(val, skips, depth))];
				end
			end
			value = strcat([value, empty]);
		end
		
%___________________________________________________________________________________________________ cellToStr
% Converts a cell array to string.
% * inCell	the input cell array														Cell
% # value	the string resulting from conversion of the cell array.						Str
		function value = cellToStr(inCell)
			cellSize        = size(inCell);
            value           = '';
            
            for i=1:cellSize(1)
                value = [value, ' {'];
                for j=1:cellSize(2)
                    if (j > 1),     value = [value, ' | ', ImogenRecord.valueToString(inCell{i,j})];
                    else            value = [value, ImogenRecord.valueToString(inCell{i,j})];
                    end
                end
                value = strcat(value, '}');
            end
		end

%___________________________________________________________________________________________________ structToStr
% Converts a structure to a string.
% * inStruct	the input structure.													Struct
% # value		the string resulting from conversion of the structure.					Str
		function value = structToStr(inStruct)			
			inCell = struct2cell(inStruct);
			for i=1:length(inCell)
				inCell{i} = ImogenRecord.valueToString(inCell{i});
			end

			inNames = fieldnames(inStruct);
			value = '';
			for i=1:length(inNames)
				value = [value, sprintf('< %s : %s >', inNames{i}, inCell{i}), blanks(2)];
			end
		end

%___________________________________________________________________________________________________ logicalToStr
% Converts a boolean to a string
% * inBool		the input logical value.												Logical
% # value		the string resulting from conversion of the boolean						Str
		function value = logicalToStr(inBool)
			value = ImogenRecord.numberToStr(inBool);
			strrep(value,'1','true');
			strrep(value,'0','false');
		end
			
	end%STATIC
	
end%CLASS


		