function updateUI(run)
%   This routine is responsible for updating the user as the code progresses. It returns a boolean
%   status variable that is used in the Imogen main iteration loop to signal the resultsManager
%   routine to save slice data at critical update iterations or time values.
%
%>< run		data manager object														ImogenManager
    
    time = run.time;
    %-----------------------------------------------------------------------------------------------
    % Clock first loop to and then estimate total run time
    %-----------------------------------------------------
	if time.iteration < 3
		switch time.iteration 
			case 1;	%Activate clock timer for the first loop	
				tic; 
			case 2;	%Stop clock timer and use the elapsed time to predict total run time
				elapsedMins = toc / 60.0;
				run.save.logPrint('\nFirst loop completed successfully in %0.5g minutes.', ...
                                  elapsedMins);
				
				if (time.iterPercent > time.timePercent)
					timeRemaining = elapsedMins*time.ITERMAX;
				else
					timeRemaining = elapsedMins*ceil(time.TIMEMAX/time.time);
				end

				dTimeNum    = timeRemaining / 1440.0;
				finTime     = now + dTimeNum;
				expComplete = datestr( finTime , 'HH:MM PM');
                
				if ( floor(finTime) - floor(now) >= 1.0 )
					expComplete = [expComplete ' on ' datestr( finTime, 'mmm-dd')];
				end
				run.save.logPrint('\n\tExpected completion time: %s', expComplete);
				
				dDays       = floor(timeRemaining/1440);
				dHours      = floor( (timeRemaining - dDays*1440)/60 );
				dMinutes    = ceil( (timeRemaining - dDays*1440 - dHours*60) );
				run.save.logPrint('\n\tWith a run time of: [%g days | %g hours | %g minutes]\n', ...
                                  dDays, dHours, dMinutes);
		end
	end

    run.save.updateDataSaves();

    %-----------------------------------------------------------------------------------------------
    % Update UI for critical loops
    %-----------------------------        
	if (run.save.updateUI)
        
		if (timePercent > iterPercent)
            infoStr = sprintf('Time: %0.5g of %0.5g', time.time, time.TIMEMAX);
            compPer = run.timePercent;
		else
			infoStr = sprintf('Iteration: %0.5g of %0.5g', time.iteration, time.ITERMAX);
			compPer = run.iterPercent;
		end
        
        %--- Prepare and display the UI update string ---%
        cTime   = now;
        curTime = strcat(datestr( cTime , 'HH:MM PM'),' on',datestr( cTime, ' mm-dd-yy'));
        run.save.logPrint('[[ %0.3g%% | %s |  %s ]]\n', compPer, infoStr, curTime);
	end
	
	run.abortCheck();
end