classdef ImageManager < handle
    % The manager class responsible for handling image related actions (primarily saving). This is a
    % singleton class to be accessed using the getInstance() method and not instantiated directly.
    
    
    %===================================================================================================
    properties (Constant = true, Transient = true) %							C O N S T A N T	 [P]
        FIELDS = {'mass','ener','momX','momY','momZ','magX', ... % The possible image types
            'magY','magZ','grav','spen','pGas', 'pTot',...
            'mach','speed','velX','velY','velZ'};
    end%CONSTANT
    
    %===================================================================================================
    properties (SetAccess = public, GetAccess = public, Transient = true) %			P U B L I C  [P]
        ACTIVE;		% Specifies if image saving is active for a run.				boolean
        INTERVAL;	% Interval between image saves.									int
        COLORMAP;	% Colormap to use for the images.								double
        frame;		% Frame index number for the next image save operation.			int
        
        
        %--- Specify whether or not the saving of an image type is allowed ---%		boolean
        mass;   ener;   spen;
        momX;   momY;   momZ;
        magX;   magY;   magZ;
        pGas;   pTot;	grav;
        velX;   velY;   velZ;
        mach;   speed;
        
        logarithmic;  % Contains fields to save as logarithmic images.              struct
        
    end%PUBLIC
    
    %===================================================================================================
    properties (SetAccess = public, GetAccess = protected) %				   P R O T E C T E D [P]
        pColordepth;	% # of color values to use in creating the colormaps		int
        parent;			% parent manager											ImogenArray
        pActive;		% active image saving slices
    end %PROTECTED
    
    
    
    
    
    %===================================================================================================
    methods %																	  G E T / S E T  [M]
    end%GET/SET
    
    %===================================================================================================
    methods (Access = public) %														P U B L I C  [M]
        
        %___________________________________________________________________________________________________ preliminary
        % Preliminary actions setting up image saves for the run. Determines which image slices to save.
        function preliminary(obj)
            obj.pActive = obj.parent.save.ACTIVE(4:6);
            if ~any(obj.pActive)
                [minval, mindex] = min(obj.parent.gridSize); %#ok<ASGLU>
                switch mindex
                    case 1;		obj.pActive(3) = true;
                    case 2;		obj.pActive(2) = true;
                    case 3;		obj.pActive(1) = true;
                end
            end
        end
        
        %___________________________________________________________________________________________________ activate
        % Activates the ImageManager by determining if any images have been enabled for saving.
        function activate(obj)
            obj.ACTIVE = false;
            for i=1:length(obj.FIELDS),	obj.ACTIVE = ( obj.ACTIVE || obj.(obj.FIELDS{i}) ); end
        end
        
        %___________________________________________________________________________________________________ getColormap
        function createColormap(obj, type, colordepth)
            obj.pColordepth = colordepth;
            switch (type)
                case 'jet';     obj.COLORMAP = jet(colordepth);
                case 'hot';     obj.COLORMAP = hot(colordepth);
                case 'bone';    obj.COLORMAP = bone(colordepth);
                case 'copper';  obj.COLORMAP = copper(colordepth);
                otherwise;      obj.COLORMAP = jet(colordepth);
            end
        end
        
        %___________________________________________________________________________________________________ imageSaveHandler
        % Handles saving of images to files.
        function imageSaveHandler(obj, mass, mom, ener, mag, grav)
            if ~( obj.ACTIVE && ~mod(obj.parent.time.iteration, obj.INTERVAL) ); return; end
            for i=4:6 % For each possible 2D slice
                if ~obj.pActive(i-3), continue; end
                
                for j=1:length(ImageManager.FIELDS)
                    if ~obj.(ImageManager.FIELDS{j}); continue; end
                    
                    switch ImageManager.FIELDS{j}
                        case 'mass'
                            array = obj.parent.save.getSaveSlice(mass.array,i);
                        case 'momX'
                            array = obj.parent.save.getSaveSlice(mom(1).array,i);
                        case 'momY'
                            array = obj.parent.save.getSaveSlice(mom(2).array,i);
                        case 'momZ'
                            array = obj.parent.save.getSaveSlice(mom(3).array,i);
                        case 'ener'
                            array = obj.parent.save.getSaveSlice(ener.array,i);
                        case 'magX'
                            array = obj.parent.save.getSaveSlice(mag(1).array,i);
                        case 'magY'
                            array = obj.parent.save.getSaveSlice(mag(2).array,i);
                        case 'magZ'
                            array = obj.parent.save.getSaveSlice(mag(3).array,i);
                        case 'grav'
                            array = obj.parent.save.getSaveSlice(grav.array,i);
                            
                        case 'spen'
                            array = ener.array ./ mass.array;
                            array = obj.parent.save.getSaveSlice(array,i);
                            
                        case 'pTot'
                            array = pressure('total', obj.parent, mass, mom, ener, mag);
                            array = obj.parent.save.getSaveSlice(array,i);
                            
                        case 'pGas'
                            array = pressure('gas', obj.parent, mass, mom, ener, mag);
                            array = obj.parent.save.getSaveSlice(array,i);
                            
                        case 'mach'
                            array = getMach(mass, mom, ener, mag, obj.parent.GAMMA);
                            array = obj.parent.save.getSaveSlice(array,i);
                            
                        case 'speed'
                            array = sqrt(getVelocitySquared(mass,mom));
                            array = obj.parent.save.getSaveSlice(array,i);
                            
                        case 'velX'
                            array = obj.parent.save.getSaveSlice(mom(1).array./mass.array,i);
                            
                        case 'velY'
                            array = obj.parent.save.getSaveSlice(mom(2).array./mass.array,i);
                            
                        case 'velZ'
                            array = obj.parent.save.getSaveSlice(mom(3).array./mass.array,i);
                            
                    end
                    
                    name = ImageManager.FIELDS{j};
                    
                    if isfield(obj.logarithmic,ImageManager.FIELDS{j})
                        logArray = ImageManager.findLog(array); % Convert array to natural log.
                        logArray = obj.parent.save.getSaveSlice(logArray,i);
                        logName = sprintf('log_%s',ImageManager.FIELDS{j});                                  
                        if labindex == 1
                            obj.saveImage(logArray, logName, obj.parent.save.SLICELABELS{i});                
                        end                       
                    end
                    
                    if labindex == 1
                        obj.saveImage(array, name, obj.parent.save.SLICELABELS{i});                
                    end
                    labBarrier(); % Block until lab 1 finishes writing images
                end
            end
            obj.frame = obj.frame + 1;
        end

        
        
    end%PUBLIC
    
    %===================================================================================================
    methods (Access = protected) %											P R O T E C T E D    [M]
        
        %___________________________________________________________________________________________________ saveImage
        %   This helper routine is responsible for saving image files. The default format here is an 8bit
        %   RGB png file.
        % array         array slice to write to an image file                           double  [nx ny]
        % name			name of the array for folder and file                           str     *
        % sliceType     slice type identifier to include in file name (eg. XY or YZ)    str     *
        function saveImage(obj, array, name, sliceType)
            minVal = min(min(array));
            maxVal = max(max(array));
            rescaledArray = obj.pColordepth * (array' - minVal) / (maxVal - minVal);
            
            iterStr = obj.parent.paths.iterationToString(obj.frame);
            fileName = strcat(name,'_',sliceType,'_',iterStr,'.png');
            filePathName = strcat(obj.parent.paths.image,filesep,name,filesep,fileName);
            imwrite(rescaledArray, obj.COLORMAP, filePathName, 'png', 'CreationTime', ...
                datestr(clock,'HH:MM mm-dd-yyyy'));
        end
        
        %___________________________________________________________________________________________________ ImageManager
        % Creates a new ImageManager instance.
        function obj = ImageManager()
            obj.frame = 0;
            for i=1:length(obj.FIELDS),     obj.(obj.FIELDS{i}) = false;    end
        end
        
    end%PROTECTED
    
    %===================================================================================================
    methods (Static = true) %													  S T A T I C    [M]
        
        %___________________________________________________________________________________________________
        % Find natural log of absolute value, replacing infinities with next greater finite minimum.
        function result = findLog(array)           
            newMin = min(log(abs(nonzeros(array))));
            result = log(abs(array));    
            infinities = isinf(result);
            [I,J] = find(infinities == true);
            for k = 1:length(I);
                result(I(k),J(k)) = newMin;
            end      
        end        
        
        
        %___________________________________________________________________________________________________ getInstance
        % Accesses the singleton instance of the ImageManager class, or creates one if none have
        % been initialized yet.
        function singleObj = getInstance()
            persistent instance;
            if isempty(instance) || ~isvalid(instance)
                instance = ImageManager();
            end
            singleObj = instance;
        end
        
    end%STATIC
end
