classdef Analyzer < handle
    % Base class for experiment analysis programs. 

    %===================================================================================================
    properties (Access = private) %                                                    P R I V A T E [P]

        pData;      % instance of the DataHandler class                             DataHandler
        pWindow;    % instance of the MultiPlot class                               MultiPlot
        varList;    % Cell of strings, input as a constructor argument.             cell        

    end %PROTECTED    

    %===================================================================================================   
    properties (Dependent)%= true, SetAccess = private) %                          D E P E N D E N T [D]

        path;       % Path to location of data files.                               str
        data;       % Storage for current data structure.                           struct
        matNames    % Data structure names, stored as strings.                      cell
        index;      % Index of current data structure.                              double
        fields;     % Use to directly access struct specified by fieldNames         struct
        fieldNames; % Names (as strings) of obj.vars fields to be analyzed.         cell
        fieldIndex; % Current index of obj.fields.                                  double
        dimensions; % Number of dimensions of simulation data being analyzed.       str
        preCalc;    % Extra calculation before plotting. Ex:'array = log(abs(x))'.  str
        plotType;   % Plotting type enumeration                                     str 
        saveImages; % True if image saving desired.                                 bool    

    end %DEPENDENT

    %===================================================================================================
    methods %                                                                         G E T / S E T [M]

    %___________________________________________________________________________________________________ GS path                
        function result = get.path(obj);    result = obj.pData.path;    end
        function set.path(obj,value);       obj.pData.path = value;     end        
    %___________________________________________________________________________________________________ GS data
        function result = get.data(obj);    result = obj.pData.vars;    end
        function set.data(obj,value);       obj.pData.vars = value;     end
    %___________________________________________________________________________________________________ GS fieldNames                        
        function result = get.matNames(obj); result = obj.pData.matNames;   end
        function set.matNames(obj,value);     obj.pData.matNames = value;   end
    %___________________________________________________________________________________________________ GS index
        function result = get.index(obj);   result = obj.pData.N1;      end
        function set.index(obj,value);      obj.pData.N1 = value;       end        
    %___________________________________________________________________________________________________ GS fields                        
        function result = get.fields(obj); result = obj.pData.fields;   end
    %___________________________________________________________________________________________________ GS fieldNames                        
        function result = get.fieldNames(obj); result = obj.pData.fieldNames;   end
        function set.fieldNames(obj,value);     obj.pData.fieldNames = value;   end
    %___________________________________________________________________________________________________ GS index
        function result = get.fieldIndex(obj);   result = obj.pData.N2;      end
        function set.fieldIndex(obj,value);      obj.pData.N2 = value;       end          
    %___________________________________________________________________________________________________ GS preCalc                    
        function result = get.dimensions(obj); result = obj.pData.dimensions;   end
        function set.dimensions(obj,value);    obj.pData.dimensions = value;    end
    %___________________________________________________________________________________________________ GS preCalc                    
        function result = get.preCalc(obj); result = obj.pData.preCalc; end
        function set.preCalc(obj,value);    obj.pData.preCalc = value;  end
    %___________________________________________________________________________________________________ GS plotType                   
        function result = get.plotType(obj);    result = obj.pWindow.plotType;  end
        function set.plotType(obj,value);   obj.pWindow.plotType = value;       end        
    %___________________________________________________________________________________________________ GS saveImages                        
        function result = get.saveImages(obj); result = obj.pWindow.saveImages; end
        function set.saveImages(obj,value); obj.pWindow.saveImages = value;     end

    end %GET/SET
    
    %===================================================================================================
    methods (Access = public) %                                                         P U B L I C  [M]
       
    %___________________________________________________________________________________________________ Analyzer
    % Constructor method, creates required instances of data and windowing objects.           
        function obj = Analyzer(varargin)    

            % Create required objects
            obj.pWindow = MultiPlot();
            obj.pData   = DataHandler();
            
            % Modify object properties
            if ~isempty(varargin);   obj.pData.matNames = varargin;   end            
            obj.plotType   = MultiPlot.PLOT_SURF2D;
            obj.saveImages = false;
                 
        end          

    %___________________________________________________________________________________________________ runAnalyzer
    % Send data and plotting objects to plotting method.        
        function runAnalyzer(obj)   
            obj.pWindow.plotFigures(obj.pData);  
        end
                
    end%PUBLIC

    %===================================================================================================
    methods (Static = true) %                                                         S T A T I C    [M] 

    %___________________________________________________________________________________________________ delete
    % Remove temporary objects when oject is destroyed.
        function delete()
            clear;
        end
            
    end % STATIC
    
end %CLASS
