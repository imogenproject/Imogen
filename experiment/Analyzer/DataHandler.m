 classdef DataHandler < handle
% Manage data to be analyzed via Analyzer class.
	
%===================================================================================================
    properties (SetAccess = public, GetAccess = public) %							P U B L I C  [P]
        vars;              % struct     Actual data structure at current index.   
        fields;            % array      Actual data field at current index.
        preCalc;           % str		Additional calculation to perform before plotting.
        fieldNames;        % cell		Names (as strings) of obj.vars fields to be analyzed.
        matNames;          % cell		Data structure names, stored as strings.
        N1;                % double     Temporary storage for obj.matNames index count.
        N2;                % double     Temporary storage for obj.fieldNames index count.
        path;              % str        Path to actual data files.
        dimensions = '2D'; % str		Specifies data to load with storeMatNames.
    end %PUBLIC
%===================================================================================================
    properties (SetAccess = protected, GetAccess = protected) %				   P R O T E C T E D [P]        
        pOriginalDirectory;    % char		Beginning location.
        pDataDirectory;        % char		Location of .mat files to be analyzed.
    end %PROTECTED
		
%===================================================================================================
    methods %																	  G E T / S E T  [M]

%___________________________________________________________________________________________________ GS vars
% Associate obj.matNames index with actual data, and execute any other preparatory calculations.
        function result = get.vars(obj)
            
            % Navigate to directory if needed, actuate obj.vars.
            if ~isempty(obj.pDataDirectory)
                cd(obj.pDataDirectory);
                load(obj.matNames{obj.N1});
                structName = who('sx*');
                result = eval(structName{:});
                cd(obj.pOriginalDirectory);
            else
                result = evalin('base',obj.matNames{obj.N1});
            end       
        
        end
%___________________________________________________________________________________________________ GS fields
% Associate obj.fieldNames index with actual data, and execute any other preparatory calculations.        
        function array = get.fields(obj)
            
            if isempty(obj.fieldNames{obj.N2}), array = []; return; end
            
            % Calculate magnitude of velocity.
            if strcmpi(obj.fieldNames{obj.N2},'speed') == true;
                array = sqrt(obj.vars.momX .^2 + obj.vars.momY .^2 ...
                                + obj.vars.momZ .^2) ./ obj.vars.mass;

            % Calculate magnitude of momentum.
            elseif strcmpi(obj.fieldNames{obj.N2},'magMom') == true;
                array = sqrt(obj.vars.momX .^2 + obj.vars.momY .^2 ...
                                                 + obj.vars.momZ .^2);
                                             
            % Calculate magnitude of magnetic field.
            elseif strcmpi(obj.fieldNames{obj.N2},'magMF') == true;
                array = sqrt(obj.vars.magX .^2 + obj.vars.magY .^2 ...
                                                + obj.vars.magZ .^2);            
            % Or just return array.
            else
                array = eval(['obj.vars.' obj.fieldNames{obj.N2}]);
            end
                        
            % Perform additional operations before returning data.
            if ~isempty(obj.preCalc)
                x = array; %#ok<NASGU>
                eval([obj.preCalc ';']);
            end
            
        end
%___________________________________________________________________________________________________ GS N1
% Determine obj.matNames indices, set to 1 if empty.        
        function result = get.N1(obj)
                
            if isempty(obj.N1)
                result = 1;
            else
                result = obj.N1;         
            end
        end
%___________________________________________________________________________________________________ GS N2
% Determine obj.fieldNames indices, set to 1 if empty.               
        function result = get.N2(obj)  
            if isempty(obj.N2)
                result = 1; 
            else
                result = obj.N2;
            end       
        end              
%___________________________________________________________________________________________________ GS fieldNames
% Determine obj.fieldNames, set to {''} if empty.       
        function result = get.fieldNames(obj)           
            if isempty(obj.fieldNames)
                result = {''};
            else
                result = obj.fieldNames;
            end
        end       
%___________________________________________________________________________________________________ GS path
% Determine path, set to [] if empty.    
        function set.path(obj,value)         
            if isempty(value)
                obj.path = [];
            elseif strcmpi(value,'workspace')           
                importWorkspace(obj);
                obj.path = value;
            else
                storeMatNames(obj,value);
                obj.path = value;
            end                
        end
               
    end%GET/SET
	
%===================================================================================================	
	methods (Access = protected) %											P R O T E C T E D    [M]

%___________________________________________________________________________________________________ importWorkspace
% Import and sort variable names from the base workspace.    
        function  importWorkspace(obj)
            
            buffer = evalin('base','who(''sx_*'')');
            buffer = sort(buffer);
            vlen = length(buffer);
            index = 1;
            obj.matNames = cell(1,vlen);
            
            for i=1:vlen              
                startRes = strfind(buffer{i},'START');                
                if ~isempty(startRes)
                    obj.matNames{1} = buffer{i};
                    index = 2;
                end          
                finalRes = strfind(buffer{i},'FINAL');              
                if ~isempty(finalRes)
                    obj.matNames{end} = buffer{i};
                end              
            end
            
            for i=1:vlen        
                startRes = strfind(buffer{i},'START');
                finalRes = strfind(buffer{i},'FINAL');              
                if isempty(startRes) && isempty(finalRes)
                    obj.matNames{index} = buffer{i};
                    index = index + 1;
                end                
            end
            
        end
        
%___________________________________________________________________________________________________ storeMatNames
% Navigate to folder specified in obj.path, sort and import variable names into obj.matNames, 
% then return to original directory.       
        function storeMatNames(obj,path)     

            % Navigate to directory.
            obj.pOriginalDirectory = pwd;
            cdPath = sprintf('cd %s',path);
            eval(cdPath);

            % Gather and sort .mat file names.
            matfiles = dir([obj.dimensions '*.mat']);
            vlen = length(matfiles);
            index = 1;
            obj.matNames = cell(1,vlen);

            for i=1:vlen
                startRes = strfind(matfiles(i).name,'START');
                if ~isempty(startRes)
                    obj.matNames{1} = matfiles(i).name;
                    index = 2;
                end
                finalRes = strfind(matfiles(i).name,'FINAL');
                if ~isempty(finalRes)
                    obj.matNames{end} = matfiles(i).name;
                end
            end

            for i=1:vlen
                startRes = strfind(matfiles(i).name,'START');
                finalRes = strfind(matfiles(i).name,'FINAL');
                if isempty(startRes) && isempty(finalRes)
                    obj.matNames{index} = matfiles(i).name;
                    index = index + 1;
                end
            end

            % Return to original directory.
            obj.pDataDirectory = pwd;
            cd (obj.pOriginalDirectory);

        end          
       
    end%PROTECTED	
end%CLASS