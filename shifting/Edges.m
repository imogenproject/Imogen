classdef Edges < handle
% Class to handle edge related settings and functionality.
	
%===================================================================================================
	properties (Constant = true, Transient = true) %							C O N S T A N T	 [P]
		DIMENSION = {'x','y','z'};
		FIELDS = {'lower','upper'};
    end%CONSTANT
	
%===================================================================================================
    properties (SetAccess = public, GetAccess = public) %							P U B L I C  [P]
		ACTIVE;
		TOLERANCE;
		lower;
		upper;
		
    end %PUBLIC

%===================================================================================================
    properties (SetAccess = protected, GetAccess = protected) %				   P R O T E C T E D [P]
		pIndex;
		pSlices;
    end %PROTECTED
	
	
	
	
	
%===================================================================================================
    methods %																	  G E T / S E T  [M]
	end%GET/SET
	
%===================================================================================================
    methods (Access = public) %														P U B L I C  [M]
		
%___________________________________________________________________________________________________ Edges
		function obj = Edges(bcModes, array, tolerance)
			if isempty(array)
				error(['Edges:EdgeError: Edge conditions set on empty array. ', ... 
						'Missing data may cause shifting errors.']);
			end
			
            %--- Initialize edge parameters ---%
			N               = size(array);
			dim             = ndims(array);
			obj.TOLERANCE   = tolerance;
			
			if (dim > 2)
				obj.pSlices = ones(3,2);
				obj.pIndex  = {1:N(1),1:N(2),1:N(3)};
			else
				obj.pSlices = ones(2,2);
				obj.pIndex  = {1:N(1),1:N(2)};
			end
			obj.pSlices(:,2) = size(array);
			
            %--- Store Edge ICs ---%
			obj.ACTIVE      = false(2, 3);
            for n=1:2
                for i=1:dim
                    iIndex      = obj.pIndex;
                    switch bcModes{n, i}

                        case ENUM.BCMODE_TRANSPARENT
                            obj.ACTIVE(n,i) = true;

                            if (n == 1)
                                iIndex{i}   = 1;
                                field       = 'lower';
                            else
                                iIndex{i}   = size(array,i);
                                field       = 'upper';
                            end
                            
                            obj.(field).(Edges.DIMENSION{i}) = array(iIndex{:});
                            
                            if ~isa(array,'double') %r2009b: iscodistributed
                                obj.(field).(Edges.DIMENSION{i}) ...
                                                         = gather(obj.(field).(Edges.DIMENSION{i}));
                            end

                            obj.(field).(Edges.DIMENSION{i}) = squeeze(obj.(field).(Edges.DIMENSION{i}));

                        case ENUM.BCMODE_FADE
                            obj.ACTIVE(n,i) = true;

                            if (n == 1)
                                iIndex{i}   = 1:20;
                                field       = 'lower';
                            else
                                iIndex{i}   = (N(i)-19):N(i);
                                field       = 'upper';
                            end
                            
                            obj.(field).(Edges.DIMENSION{i}) = array(iIndex{:});
                            
                        case ENUM.BCMODE_WALL
                            obj.ACTIVE(n,i) = true;
                            
                            if (n == 1)
                                iIndex{i} = 1;
                                field     = 'lower';
                            else
                                iIndex{i} = N(i);
                                field     = 'upper';
                            end
                            
                            obj.(field).(Edges.DIMENSION{i}) = array(iIndex{:});
                    end
                end
            end
		end

%___________________________________________________________________________________________________ getEdge
		function result = getEdge(obj, upper, dim, array, bcType)
			
			upper = upper + 1;
			
			switch bcType
				
				case ENUM.BCMODE_TRANSPARENT
					iIndex		= obj.pIndex;
					iIndex{dim} = obj.pSlices(dim,upper);

					newEdge     = array(iIndex{:});
                    if ~isa(array,'double') %r2009b: iscodistributed
                        newEdge = gather(newEdge);	
                    end 
					newEdge     = squeeze(newEdge); 
					oldEdge     = obj.(Edges.FIELDS{upper}).(Edges.DIMENSION{dim});

					delta       = min(abs(newEdge - oldEdge),obj.TOLERANCE);
					signTest    = (newEdge - oldEdge) > 0;
					result      = oldEdge + (signTest - ~signTest) .* delta;

					obj.(Edges.FIELDS{upper}).(Edges.DIMENSION{dim}) = result; %Update edge
					
                
				case {ENUM.BCMODE_FADE, ENUM.BCMODE_WALL}
					result = obj.(Edges.FIELDS{upper}).(Edges.DIMENSION{dim});
			end
		end
		
	end%PUBLIC
	
%===================================================================================================	
	methods (Access = protected) %											P R O T E C T E D    [M]
	end%PROTECTED
	
end%CLASS
