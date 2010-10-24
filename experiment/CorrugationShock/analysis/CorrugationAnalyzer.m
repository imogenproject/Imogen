classdef CorrugationAnalyzer < Analyzer
% Analyzer class for corrugation shock problems.
	
%===================================================================================================
	properties (Constant = true, Transient = true) %							C O N S T A N T	 [P]        
		VERSION = 1.3;
    end%CONSTANT
	
%===================================================================================================
    properties (Dependent = true) %											  D E P E N D E N T  [P]         
        fftData;
        growthData;
        crossingTimes;        
    end %DEPENDENT

%===================================================================================================
    properties (SetAccess = protected, GetAccess = protected) %				   P R O T E C T E D [P]
       pfftData;            % stores fft data for each slice                            double(?,?)
       pDataIndex;          % next index into the data arrays                           double
       pIterations;         % iteration of each data slice                              double(?)
       pSimulationTimes;    % simulated time at each data slice                         double(?)
       pGrowthData;         % growth of the instability using Stone & Edelman format    double(?)        
       pTimesteps;          % history of timesteps for the run                          double(?)
       pTreadmills;         % treadmill activity for the run                            double(?)
       pGridProfile;        % profile of dGrid along the x direction                    double(?)
       pMagnetAmpSq0;       % initial pre-shock magnetic field amplitude squared        double
       pAlfvenVelocity;     % pre-shock Alfven velocity at earliest run state           double
    end %PROTECTED
	
%===================================================================================================
    methods %																	  G E T / S E T  [M]       

%___________________________________________________________________________________________________ GS: fftData
        function result = get.fftData(obj)
            try
                result = obj.pfftData(1:obj.pDataIndex, :);
            catch MERR %#ok<*NASGU>
                result = [];
            end
        end

%___________________________________________________________________________________________________ GS: growthData
        function result = get.growthData(obj)
            result = obj.pGrowthData;
        end
        
%___________________________________________________________________________________________________ GS: crossingTimes
        function result = get.crossingTimes(obj)
            %--- Determine crossing times based on stored simulation times ---%
            %       Here I define the crossing distance as is it is defined in the Stone & Edeleman 
            %       paper by the parameter L = 0.01, which is also used to set the grid size in the 
            %       CorrugationShockInitializer.
            if isempty(obj.pAlfvenVelocity)
                result = [];
            else
                crossingDistance = 0.01; %See Stone & Edelman
                result = obj.pSimulationTimes ./ (crossingDistance/obj.pAlfvenVelocity);
            end
        end   
    
    end%GET/SET
	
%===================================================================================================
    methods (Access = public) %														P U B L I C  [M]
        
%___________________________________________________________________________________________________ CorrugationAnalyzer
% Initializes properties, overloads desired properties of Analyzer class.        
        function obj = CorrugationAnalyzer(filename)
            obj.pDataIndex = 1;
            if (nargin > 0 && ~isempty(filename));
            	obj.fromFile(filename);
            else
            	obj.dimensions = '3D';
            	obj.path = 'workspace'; % Options: 'workspace' or '/path/to/matfiles/'.
                obj.plotType = 'none';    	
            end
        end
        
%___________________________________________________________________________________________________ run
% Perform all desired functions in one step.        
        function run(obj)
            analyzeFiles(obj);
            plotResults(obj);
        end       

%___________________________________________________________________________________________________ analyzeFiles
        function analyzeFiles(obj)

            %--- Add sequential files ---%
            for i = 1:length(obj.matNames)
                obj.index = i;
                fprintf('Analyzing %s...........',obj.matNames{i});
                obj.addAnalysisPoint(obj.data);
                fprintf('complete\n');
            end 
            fprintf('Analysis complete\n\n');
        end
        
%___________________________________________________________________________________________________ addAnalysisPoint
        function addAnalysisPoint(obj, inputData)
                
            %--- Initialization ---%
            dataAsArrays = resultsToData(inputData);
            
            %--- Store valuable information ---%
            if isempty(obj.pGridProfile)
                if numel(inputData.dGrid{1}) > 1
                    obj.pGridProfile = inputData.dGrid{1}(:,1,1);
                else
                    obj.pGridProfile = inputData.dGrid{1}*ones(size(inputData.mass,1),1);
                end
            end
            
            %-- FINAL preferred data ---%
            if isempty(obj.pIterations) || inputData.time.iteration > max(obj.pIterations)
               obj.pTreadmills = inputData.tread.history;
               obj.pTimesteps  = inputData.time.history;
            end
            
            %--- START preferred data ---%
            if isempty(obj.pIterations) || inputData.time.iteration < min(obj.pIterations)
                obj.pMagnetAmpSq0 = inputData.magX(1,1,1)^2 + inputData.magY(1,1,1)^2 ...
                                        + inputData.magZ(1,1,1)^2;
                obj.pAlfvenVelocity   = maxFinderND( sqrt(getMagneticSquared(dataAsArrays.mag) ...
                                        ./ inputData.mass) );
            end
            
            obj.pSimulationTimes(obj.pDataIndex)  = sum(inputData.time.history);
            obj.pIterations(obj.pDataIndex)       = inputData.time.iteration;
            
            
            %--- Analysis ---%
            obj.fftAnalysis(inputData);
            obj.growthAnalysis(dataAsArrays);
            
            obj.pDataIndex = obj.pDataIndex + 1; % Step data index for next input.
        end 
    
%___________________________________________________________________________________________________ plotResults
        function plotResults(obj, visible, print)    
            
            %--- Initialization ---%
            if nargin < 3 || isempty(print);    print = false; end
            
            if nargin < 2 || isempty(visible);  visible = 'on'; 
            else
                if (visible); visible = 'on'; else visible = 'off'; end
            end
            
            width = 2;
            
            %---------------------------------------------------------------------------------------
            % Time analysis plot
            %-------------------
            hFig                = figure('Color',[1 1 1],'Visible',visible);
            numTicks            = min(5, length(obj.pIterations));
            crossingIndeces     = round(linspace(1,length(obj.pIterations),numTicks));
            crossingTicksVec    = ones(1,numTicks);
            crossingTicks       = cell(1,numTicks);
            for i=1:length(crossingIndeces)
                crossingTicksVec(i) = obj.pIterations(crossingIndeces(i));
                crossingTicks{i}    = num2str(obj.crossingTimes(crossingIndeces(i)),'%0.3g');
            end 
            
            %--- Time elapsed ---%
            hAxis = subplot(2,2,1);
            plot(cumsum(obj.pTimesteps),'LineWidth',width,'Color',[0 0 1]);
            set(hAxis, 'XTick', crossingTicksVec, 'XTickLabel', crossingTicks);
            xlabel('Crossing times $\left(\frac{t}{t_A}\right)$','Interpreter','latex');
            title('Elapsed Time');
            grid on;

            %--- Timestep ---%
            hAxis = subplot(2,2,2);
            plot(obj.pTimesteps,'LineWidth',width,'Color',[0 0 1]);
            set(hAxis, 'XTick', crossingTicksVec, 'XTickLabel', crossingTicks);
            xlabel('Crossing times $\left(\frac{t}{t_A}\right)$','Interpreter','latex');
            title('Timesteps');
            grid on;

            %--- Treadmill ---%
            hAxis = subplot(2,2,3);
            plot(abs(obj.pTreadmills(1:max(obj.pIterations))),'LineWidth',width,'Color',[1 0 0]);
            set(hAxis, 'XTick', crossingTicksVec, 'XTickLabel', crossingTicks);
            xlabel('Crossing times $\left(\frac{t}{t_A}\right)$','Interpreter','latex');
            title('Treadmill Activity');
            grid on;

            %--- dGrid ---%
            hAxis = subplot(2,2,4);
            plot(obj.pGridProfile,'LineWidth',width,'Color',[0 0.5 0]);
            xlabel('Cells (shock normal direction)');
            ylabel('Cell length');
            title('Grid Spacing Profile (Shock Direction)');
            grid on;
            
            if print;   printResultsFigure('time_results', hFig, {'png'}); end
            
            %---------------------------------------------------------------------------------------
            % Instability plot
            %-----------------
            hFig  = figure('Color', [1 1 1],'Visible',visible);
            
            %--- Growth parameter plot ---%
            hAxis = subplot(1,2,1); hold all;
            plot(obj.crossingTimes, obj.pGrowthData,'Color',[0,0,0],'LineWidth',width,'Marker','o');
            grid on;
            box on;
            title('Instability Growth','Interpreter','latex');
            xlabel('Crossing times $\left(\frac{t}{t_A}\right)$','Interpreter','latex');
            ylabel('$\xi$','Interpreter','latex');
            
            
            %--- FFT plot ---%
            freq = obj.pfftData(1:(obj.pDataIndex-1),:);
            hAxis = subplot(1,2,2); hold all;
            bar3(freq);
            ylabel('Crossing times $\left(\frac{t}{t_A}\right)$','Interpreter','latex');
            xlabel('Spatial frequency $\left(length^{-1}\right)$','Interpreter','latex');
            zlabel('Amplitude (normalized)','Interpreter','latex');
            grid on;
            title('FFT Amplitudes of Corrugation Over Time','Interpreter','latex');
            numFreqs = size(freq,2);
            colors   = 0.5*ones(numFreqs,3);
            colors(1:3:end,1) = 0.9;
            colors(2:3:end,1) = 0.75; colors(2:3:end,3) = 0.75;
            colors(3:3:end,3) = 0.9;
            colormap(colors);
            set(hAxis, 'XLim',[0, size(freq, 2) + 1],'Ylim',[1, size(freq, 1) + 1]);
            
            % Total length = 0.01 so min frequency 100, max frequency = length(freq)/0.01
            vals    = linspace(100, size(freq,2)/0.01, size(freq,2));
            indeces = 1:5:length(vals);
            if (indeces(end) ~= length(vals)),   indeces(end+1) = length(vals); end
           
            ticks = cell(1, length(indeces));
            for i=1:length(ticks),  ticks{i} = num2str(vals(indeces(i)),'%0.5g'); end
            set(hAxis, 'XTick', indeces, 'XTickLabel', ticks);
            set(hAxis, 'YTick', crossingIndeces, 'YTickLabel', crossingTicks);
            view(hAxis,[143 28]);
            
            if print;   printResultsFigure('instability_results', hFig, {'png'}); end
        end    
         
%___________________________________________________________________________________________________ fromFile 
        function fromFile(obj, filename) 
            res = load(filename); 
            if ~isfield(res,'corr'); 
                error('Unable to load file. Data improperly saved.'); 
            end 
             
            obj.pfftData        = res.corr.pfftData; 
            obj.pDataIndex      = res.corr.pDataIndex; 
            obj.pIterations     = res.corr.pIterations; 
            obj.pGrowthData     = res.corr.pGrowthData; 
            obj.pTimesteps      = res.corr.pTimesteps; 
            obj.pTreadmills     = res.corr.pTreadmills; 
            obj.pGridProfile    = res.corr.pGridProfile; 
            obj.pMagnetAmpSq0   = res.corr.pMagnetAmpSq0; 
            obj.pSimulationTimes = res.corr.pSimulationTimes; 
            obj.pAlfvenVelocity = res.corr.pAlfvenVelocity; 
            fprintf('Data loaded and ready.\n'); 
        end 
         
%___________________________________________________________________________________________________ 
        function toFile(obj, filename) 
            corr.pfftData       = obj.pfftData; 
            corr.pDataIndex     = obj.pDataIndex; 
            corr.pIterations    = obj.pIterations; 
            corr.pGrowthData    = obj.pGrowthData; 
            corr.pTimesteps     = obj.pTimesteps; 
            corr.pTreadmills    = obj.pTreadmills; 
            corr.pGridProfile   = obj.pGridProfile; 
            corr.pMagnetAmpSq0  = obj.pMagnetAmpSq0; 
            corr.pSimulationTimes = obj.pSimulationTimes; 
            corr.pAlfvenVelocity = obj.pAlfvenVelocity;  %#ok<STRNU>
            save(filename,'corr'); 
            fprintf('Data saved to: %s\n',filename); 
        end          
         
    end%PUBLIC
    
%===================================================================================================	
	methods (Access = protected) %											P R O T E C T E D    [M]

%___________________________________________________________________________________________________ growthAnalysis
        function growthAnalysis(obj, data)
        % Analyzes the growth of the corrugation instability
                
            %--- Get the x-compression term ---%
            velocity = data.mom(1).dataClone();
            velocity.array =  velocity.array ./ data.mass.array;
            compress = velocity.calculate2PtDerivative(1,data.run.dGrid{1});
            compress = min(compress, 0);

            %--- Calculate growth parameter ---%
            %       The growth parameter is calculated according to the Stone & Edeleman '93 paper.
            N = data.mass.gridSize;
            growthParam = compress .* data.mag(3).array .* data.mag(3).array;
            growthParam = sumND( growthParam(25:(N(1)-25),:,:) ) ...
                                                            ./ sumND( compress(25:(N(1)-25),:,:) );
            growthParam = growthParam ./ obj.pMagnetAmpSq0;
            obj.pGrowthData(obj.pDataIndex) = log(growthParam);
        end
        
%___________________________________________________________________________________________________ fftAnalysis
        function fftAnalysis(obj, inputData)
        % Calculates the FFT for the input data and adds the results to the object's data for later
        % plotting.
        
            %--- Initialization ---%
            grid            = size(inputData.mass);
            freq            = zeros(1,grid(2));
            shockSlice      = squeeze(inputData.mass(ceil(grid(1)/2),:,:));
            
            %--- Y-directed FFTs ---%
            for j=1:grid(2)
                freq = freq + abs( fft(shockSlice(j,:) - mean(shockSlice(j,:))) );
            end

            %--- Z-directed FFTs ---%
            for k=1:grid(3)
                freq = freq + abs( fft(shockSlice(:,k) - mean(shockSlice(:,k))) )';
            end

            %--- Cleanup frequency results ---%
            %       FFT is symmetric about 0 and shifted by pi/2, and in this case the mean was 
            %       subtracted off during fft analysis, so the 0 mode amplitude will always be 0. 
            %       Hence, a range of 2:(N/2+1) gives the appropriate range for display. Also, the
            %       amplitudes are renormalized so that the total spectrum amplitude sums to 1.
            maxFreq = floor(grid(2)/2) + 1;
            freq    = freq(2:maxFreq);
            freq    = freq/sum(freq);
            
            %--- Add inputData to data object ---%
            if isempty(obj.pfftData)
                obj.pfftData                        = zeros(25,length(freq));
                
            elseif (obj.pDataIndex > size(obj.pfftData,1))
                largerArray                         = zeros(size(obj.pfftData,1)+25, ...
                                                                              size(obj.pfftData,2));
                largerArray(1:(obj.pDataIndex-1),:) = obj.pfftData;
                obj.pfftData                        = largerArray;
            end
            
            obj.pfftData(obj.pDataIndex,:)   = freq;
        end    
    
    end%PROTECTED
		
%===================================================================================================	
	methods (Static = true) %													  S T A T I C    [M]
    end%STATIC

end%CLASS
