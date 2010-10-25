function plot3DArray(array,contrast,scaleMax,pTag,sliceVals)

%=== PLOT3DARRAY ===================================================================================
%   This routine plots 3D arrays by creating a 2x2 subplot grid and plotting a 3D plot as well
%   as 3 2D plots for each of the axis planes and colormap all of them according to the min and
%   max values of the array and the contrast value.
%===================================================================================================
%
%-CALLS---------------------------------------------------------------------------------------------
%   plot3DArray(array, <contrast>, <scaleMax>, <pTag>, <sliceVals>)
%-INPUTS--------------------------------------------------------------------------------------------
% array         3D or 4D input array to plot                                double  [(3) nx ny nz]
% contrast      percentage of maximum contrast                              double  #
% scaleMax      maximum color value used as: [-scaleMax,scaleMax]           double  # OR [min# max#]
% pTag          tag identifier string for the plot                          str     *
% sliceVals     array of slice locations for each axis plane                int     [sx sy sz]
%-----------------------------------------------------------------------------------------------
    
    %------------------------------------------------------------------------------------------------
    % Argument verification and parsing
    %----------------------------------
    
    if ( (nargin < 1) || isempty(array) ) %HANDLE: missing array arg
        if (nargin == 0); disp(infoStr); end
        error('No array to plot provided. Aborting attempt to plot');
    end
    
    if ( (nargin < 2) || isempty(contrast) ), contrast = 100; end %HANDLE: missing contrast arg

    if ( (nargin < 3) || isempty(scaleMax) ) %HANDLE: missing maximum scale arg
        maxA = maxFinderND(array);
        minA = minFinderND(array);
    elseif (length(scaleMax) < 2), maxA = scaleMax; minA = -1.0 * scaleMax;
    else                       minA = scaleMax(1); maxA = scaleMax(2);
    end
    scaleMax = [minA maxA];
    
    if ( (nargin < 4) || isempty(pTag) ) %HANDLE: missing plot tag arg
        for i=1:1000
            pTag = ['plot3D_' num2str(i)];
            hExists = findobj('Tag',pTag);
            if isempty(hExists); break; end
        end
    end
           
    %-----------------------------------------------------------------------------------------------
    % Verify the array is 3D
    %-----------------------
    array = squeeze(array);    
    N = size(array);
    dim = length(N); 
    if (dim < 3);       error('Array is not 3D. Unable to plot'); 
    elseif (dim == 3);  iMax = 1; nx = 2; ny = 2; names = '3D Plot';
    elseif (dim == 4);  iMax = 3; N = N(2:4); nx = 3; ny = 4; names = {'3D X Comp', '3D Y Comp', '3D Z Comp'};
    end

    if ( (nargin < 5) || isempty(sliceVals) ) %HANDLE: missing sliceVals arg
       sliceVals = floor(N / 2);
    end
    
    %-----------------------------------------------------------------------------------------------
    % Create the figure
    %------------------           
    slicePlotAbility = true; % Slice plots error out when an array max == min. Check for this.
    if( abs(maxA-minA) == 0 ) % prevent errors for case minA == maxA
        slicePlotAbility = false;
        maxA = maxA + 1e-3; 
    end 

    hFig = findobj('Tag',pTag);
    if isempty(hFig)
    hFig = figure();
        set(hFig,'Tag',pTag);
        set(hFig,'Name',pTag);
    end
    set(hFig,'Color',[1 1 1]);
    set(hFig,'UserData',array);
    
    %-----------------------------------------------------------------------------------------------
    % Fill the figure with plots
    %---------------------------
    for i=1:iMax
        
        array3D = zeros(N+1);
        % Select individual components if 4-D vector array
        if (iMax > 1) 
            array3D(1:N(1),1:N(2),1:N(3)) = squeeze(array(i,:,:,:)); offset = 4*(i-1);
            slicePlotAbility = true;
            minTest = min(array3D(:)); maxTest = max(array3D(:));
            if ( (minTest - maxTest) == 0 ); slicePlotAbility = false; end
        else array3D(1:N(1),1:N(2),1:N(3)) = array; offset = 0;
        end
        
        % --- 3D plot ---%
        hAxis = subplot(nx,ny,(offset + 1));
        hold all;
        title(names(i));
        if (slicePlotAbility == true)
            grid on;
            set(hAxis,'Box','on');
            view([-45 45]);
            xlabel('X'); ylabel('Y'); zlabel('Z');

            Np = N+1;
            xlim([1 Np(1)]); ylim([1 Np(2)]); zlim([1 Np(3)]);
            [x,y,z] = meshgrid(1:Np(1),1:Np(2),1:Np(3));
            caxis(contrast/100 * [minA maxA]);
            hPlot = slice(x,y,z,permute(array3D,[2 1 3]),sliceVals(1),sliceVals(2),sliceVals(3));
            set(hPlot,'LineStyle','none');
            caxis(contrast/100 * [minA maxA]);
        else
            text(0.1,0.5,'Single Valued Array');
        end

        % --- YZ axis slice ---%
        hAxis = subplot(nx,ny,(offset + 2));
        hold all;
        title('YZ Slice');
        xlabel('Y'); ylabel('Z');
        caxis(contrast/100 * [minA maxA]);
        hPlot = imagesc(squeeze(array3D(sliceVals(1),1:N(2),1:N(3)))');
        set(hPlot,'UserData',squeeze(array3D(sliceVals(1),1:N(2),1:N(3)))');
        set(hPlot,'ButtonDownFcn',{@plotClick_Callback,['YZ Slice at: ' mat2str(sliceVals(2:3))]});
        
        % --- XZ axis slice ---%
        hAxis = subplot(nx,ny,(offset + 3));
        hold all;
        title('XZ Slice');
        xlabel('X'); ylabel('Z');
        caxis(contrast/100 * [minA maxA]);
        hPlot = imagesc(squeeze(array3D(1:N(1),sliceVals(2),1:N(3)))');
        set(hPlot,'UserData',squeeze(array3D(1:N(1),sliceVals(2),1:N(3)))');
        set(hPlot,'ButtonDownFcn',{@plotClick_Callback,['XZ Slice at: ' mat2str([sliceVals(1) sliceVals(3)])]});
        
        % --- XY axis slice ---%
        hAxis = subplot(nx,ny,(offset + 4));
        hold all;
        title('XY Slice');
        xlabel('X'); ylabel('Y');
        caxis(contrast/100 * [minA maxA]);
        hPlot = imagesc(squeeze(array3D(1:N(1),1:N(2),sliceVals(3)))');
        set(hPlot,'UserData',squeeze(array3D(1:N(1),1:N(2),sliceVals(3)))');
        set(hPlot,'ButtonDownFcn',{@plotClick_Callback,['XY Slice at: ' mat2str(sliceVals(1:2))]});
        
    end
    
    % Show colobar for entire plot
    colorbar('peer',hAxis,[0.9164 0.1378 0.01066 0.7162]); 
    
    %--- Activate GUI ---%
    hGUI = findobj('Tag',[pTag '_GUI']);
    if isempty(hGUI)
        hGUI = Slice3D_GUI(hFig,'plot3DArray',sliceVals,array,contrast,scaleMax);
        set(hFig,'CloseRequestFcn',{@partneredCloseCB,hGUI});
    end
    
end
