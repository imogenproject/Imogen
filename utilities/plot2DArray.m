function hFig = plot2DArray(array, contrast, type)
%  This routine will plot a 2D array using a custom colormap (Min = red Mid = black MAX = white)
%  and adjust the colormap max and min values according to the percentage contrast value.
%
% array           2D input array to plot                                        double  [nx ny]
% contrast        percentage of maximum contrast                                double  #
% type            the type of 2D plot either: 'image' or 'surface'              str     *
% hFig            handle to the 2D plot figure                                  handle

    %-----------------------------------------------------------------------------------------------
    % Argument verifications and parsing
    %-----------------------------------
    if ( (nargin < 1) || isempty(array) ), error('Not enough input arguments.');  end
    if ( (nargin < 2) || isempty(contrast) ), contrast = 100; end

    if ( (nargin < 3) || isempty(type) )
        type = questdlg('Which plot type would you like?','Select Plot Type','Image','Surface','Cancel','Image');
        if strcmpi(type,'Cancel'); return; end
    end
    type = lower(type);
    
    %-----------------------------------------------------------------------------------------------
    % Verify that array is 2D
    %------------------------
    array = squeeze(array);
    N = size(array); dim = length(N); 
    vector = 1; ng = 1;
    if (dim ~= 2) 
        if ( (dim == 3) && (N(1) == 3) ) %Test to see if the array is a vector.
            vector = 3; ng = 2; N = N(2:3);
        else
            error('Array is not 2D. Unable to plot'); 
        end
    end
    
    %-----------------------------------------------------------------------------------------------
    % Create the figure and plot based on array structure
    %----------------------------------------------------
    hFig = figure;
    set(hFig,'Color',[1 1 1]);
    set(hFig,'KeyReleaseFcn',@Colormap_KeyRelease_Callback);
    for i=1:vector
        hAxis = subplot(ng,ng,i);
        if (vector > 1); resArray = squeeze(array(i,:,:));
        else resArray = array; end

        
        high = max(max(resArray)); hStr = sprintf('%0.2g',high);
        low  = min(min(resArray)); lStr = sprintf('%0.2g',low);
        if ( abs(low - high) == 0) %Prevent errors for case where high == low
            high = high + 1e-9;
        end
        
        %--- Create Colormap ---%
        midVal = (high - low)/2 + low;
        high = contrast/100*(high-midVal) + midVal;
        low = contrast/100*(low-midVal) + midVal;
        fprintf('Colormap range: [%0.2g, %0.2g].',low, high);
        caxis([low, high]);
        colormap(jet(256));
        
        
        switch type
            case 'image'
                hPlot = imagesc(resArray');
                colorbar;
                aspect = [size(resArray) 1]/min(size(resArray));
                set(hAxis,'XLim',[0.5 (N(1)+0.5)],'YLim',[0.5 (N(2)+0.5)],'YDir','normal', ...
                    'XDir','normal','PlotBoxAspectRatio',aspect,'FontSize',12);
            case 'surface'
                hPlot = surf(resArray');
                lightPlot(hFig);
                set(hAxis,'FontSize',12);
            otherwise
                warning('IMOGENUTIL:UnknownPlotType','Unknown plot type.');
        end
    end
end