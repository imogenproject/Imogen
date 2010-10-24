%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    Analyzer Class Example Run Script                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myClass = Analyzer();
    %
    % Create an instance of the Analyzer class. If you include loaded
    % structures as input arguments, i.e. Analyzer('sx_XY_01','sx_XY_02'),
    % comment out myClass.path below.
    
myClass.path = 'workspace';
    %
    % Examples:     '/Users/paulfernandez/Documents/matFiles';
    %               'workspace';
    %               
    % All .mat files in given path will be loaded. This will be overridden
    % if input arguments are provided. For Windows machines, use relative
    % path, i.e. '../../Documents/matFiles'. No spaces are allowed.
  
    
myClass.fieldNames = {'mass'};
    %
    % Example:  {'mass','speed','magMom','ener'};
    
    
myClass.dimensions = '2D';
    
    
%myClass.preCalc = 'array = DiskAnalyzer.diskUnwrap(x);';
    %
    % Example:     'array = log(abs(x))';
    % Comment out when not in use.
    
    
myClass.plotType = MultiPlot.PLOT_SURF2D;
    %
    % Examples:    @imagesc;
    %              MultiPlot.PLOT_SURF2D;
    %              MultiPlot.PLOT_IMAGE;
    
myClass.saveImages = true;
    %
    % Close figures and save as .png images.        
    
    
myClass.runAnalyzer();
    %
    % Executes class.