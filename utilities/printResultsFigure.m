function printResultsFigure(fileName, figureName, formats)
%   Prints the specified, or current, figure to a png file to a subdirectory called figures of the
%   current directory unless the current directory is itself called figures, in which case it writes
%   the file to the current directory.
%
%   fileName        name of the desired file (without extension)                        str
%   figureName      tag name for the figure to be printed                               str
%   formats         list of save formats to use when exporting. Options include:        cell
%                           'ps', 'eps', 'pdf', 'png'


    %-----------------------------------------------------------------------------------------------
    % Initialization and argument handling
    %-------------------------------------    
    if ( (nargin < 2) || (isempty(figureName)) ) %HANDLE: missing figure name
        hFig = gcf;
        disp('No figure name specified. Using current figure instead.');
    else
        if isa(figureName,'char');          hFig = findobj('Tag',figureName); 
        elseif isa(figureName,'double');    hFig = figureName;
        else
            hFig = gcf;
            disp('Invalid figure name. Printing current figure instead.');
        end 
    end
    
    if ( (nargin < 3) || isempty(formats) ), formats = {'png'}; end
    
    if ( (nargin < 1) || isempty(fileName) ) %HANDLE: missing file name
        fileName = 'results';
        disp('No file name specified. Using ''results'' instead.');
    end
    
    %-----------------------------------------------------------------------------------------------
    % Get the current folder name
    %----------------------------
    wdir = cd; % Get working directory
    test = regexp(wdir,'figures','ONCE');
    
    if isempty(test)
        if ~exist('figures','dir'); mkdir('figures'); end
        fileName = ['./figures/' fileName];
    end
    
    %-----------------------------------------------------------------------------------------------
    % Prepare and print the figure
    %-----------------------------   
    set(hFig,'Position',[0, 0, 1920, 1080]);
    
    set(hFig,'PaperType','usletter');
    set(hFig,'PaperOrientation','landscape');
    set(hFig,'PaperUnits','normalized');
    set(hFig,'PaperPosition',[0.02, 0.02, 0.98, 0.98]);
    sizeMode = '-r300';
    
    for i=1:length(formats)
        set(hFig,'PaperPositionMode','manual');
        switch (lower(formats{i}))
            case 'pdf';         print(hFig,'-dpdf',     sizeMode, fileName);
            case 'ps';          print(hFig,'-dpsc2',    sizeMode, fileName);
            case 'eps';         print(hFig,'-depsc2',   sizeMode, fileName); 
        end
    end
    
    if any(strcmpi('png',formats))
        set(hFig,'PaperPositionMode','auto');
        print(hFig,'-dpng',     sizeMode, fileName);
    end
    
end