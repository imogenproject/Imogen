classdef MultiPlot < handle
    % Loop and save figures in conjunction with Analyzer class.
    
    %===================================================================================================
    properties (Constant = true, Transient = true) %							C O N S T A N T	 [P]
        PLOT_SURF2D  = 'plot_surf2d';  % ENUMERATION: For plotting surface 2D plots.
        PLOT_IMAGE   = 'plot_image';   % ENUMERATION: For standard 2D images.  
    end%CONSTANT
    
    %===================================================================================================
    properties (SetAccess = public, GetAccess = public) %						    P U B L I C  [P]
        plotIt;                  % the handle to the plotting function             handle
        saveImages = false;      % set to true for image saving                    bool
        plotType;                % Plotting type enumeration                       str   
    end %PUBLIC

    %===================================================================================================
    methods %																	  G E T / S E T  [M]

    %___________________________________________________________________________________________________ GS plotType
    % Call Analyzer.plotSelector when plotType is set.                
        function set.plotType(obj,value) 
                plotSelector(obj,value);
                obj.plotType = value;
        end
        
    end%GET/SET
    
    %===================================================================================================
    methods (Access = public) %														P U B L I C  [M]       
        
    %___________________________________________________________________________________________________ plotSelector
    % Attaches the correct plotter function to the plotIt handle property as specified by the input type.      
        function plotSelector(obj,type)
            if ischar(type)    
                switch type
                    %-----------------------------------------------------------------------------------
                    case MultiPlot.PLOT_SURF2D
                        obj.plotIt           = @MultiPlot.createSurf2D;
                    %-----------------------------------------------------------------------------------
                    case MultiPlot.PLOT_IMAGE
                        obj.plotIt           = @MultiPlot.scMaker;
                    %-----------------------------------------------------------------------------------               
                    otherwise
                        obj.plotIt           = @MultiPlot.scMaker;
                end
            else
                obj.plotIt = type;
            end                
        end            
        
    %___________________________________________________________________________________________________ plotFigures
    % Creates a subplot for each index of data.fieldNames, looping for each index of data.MatNames               
        function plotFigures(obj,data)
            
            rect = [.1 .1 .6 .8];
            
            %--- initialize subplot limits ---%
            w = round(sqrt(length(data.fieldNames)));
            h = w;
            if (w * w) < length(data.fieldNames)
                w = w + 1;
            end
            
            %--- initialize figure window ---%
            for N1 = 1:length(data.matNames)
                data.N1 = N1;
                figure('NumberTitle','off','Units','normalized',...
                    'OuterPosition',rect);
                
                %--- create subplots ---%
                for N2 = 1:length(data.fieldNames)
                    data.N2 = N2;
                    subplot(h,w,N2);
                    
                    %--- insert plot method ---%
                    try 
                        obj.plotIt(data.fields);
                    catch %#ok<*CTCH>
                        try
                            obj.plotIt(data.vars);
                        catch
                            try
                                obj.plotIt(data);
                            catch
                                try
                                    obj.plotIt;
                                catch
                                    eval(obj.plotIt);
                                end
                            end
                        end
                    end
                    set(gca,'PlotBoxAspectRatioMode','manual');
                end
                
                name = sprintf('Iteration %d',data.vars.iter);
                set(gcf,'Name',name);
                
                if obj.saveImages   % save as .png image
                    if N1 == 1 && isdir('Analyzer_images')
                        error('Cannot overwrite existing image files: please rename or delete folder.');
                    end
                    if ~isdir('Analyzer_images')
                        mkdir Analyzer_images;
                    end
                    saveas(gca,sprintf('./Analyzer_images/%04d.png',N1),'png');
                    close;
                end
            end
        end
        
        
    end%PUBLIC
    
    %===================================================================================================
    methods (Access = protected) %											P R O T E C T E D    [M]
    end%PROTECTED
    
    %===================================================================================================
    methods (Static = true) %													  S T A T I C    [M]

    %___________________________________________________________________________________________________ createSurf2D
    % Default plotter for Analyzer class, creates a nicely rendered 2D surf.                       
         function createSurf2D(data)
            
            name = sprintf('%s at iteration %d and time %d',...
                data.fieldNames{data.N2},data.vars.iter,data.vars.time.time);
            
            %--- Create surf ---%
            surf(data.fields,'LineStyle','none');
            title (name);
            sz = size(data.fields);
            xlim ([0 sz(2)]);
            ylim ([0 sz(1)]);
            colormap(jet(256));
            view([0 0 90]) %    view([-15 75 130]);
            light('Position',[0 -2 1]);
            lightangle(-45,30);
            grid('off');
            hold('all');
            material dull;
            shading interp;
            colorbar;
            % alpha(.5);    %<--- Uncomment to set opacity:
            set(gcf,'Renderer','OpenGL');
            set(findobj(gca,'type','surface'),...
                'FaceLighting','gouraud',...
                'AmbientStrength',0.25,'DiffuseStrength',1,...
                'SpecularColorReflectance',0.25,...
                'SpecularStrength',1,'SpecularExponent',50,...
                'BackFaceLighting','lit');
         end      
        
    %___________________________________________________________________________________________________ scMaker
    % Gives imagesc a title and colorbar.                       
         function scMaker(data)
             
             imagesc(data.fields);
             title(data.fieldNames{data.N2});
             colorbar;
             
         end
         
    end%PROTECTED
    
end%CLASS