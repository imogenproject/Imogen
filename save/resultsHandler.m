function resultsHandler(run, mass, mom, ener, mag, grav)
%   Prepare and manage the storage of intermediate and then final results as a trial is run.
%
%>< run     Imogen Run manager object.                                          ImogenManager
%>< mass    Mass density object.                                                FluidArray
%>< mom     Momentum density objects.                                           FluidArray(3)
%>< ener    Ener density object.                                                FluidArray
%>< mag     Magnetic field density objects.                                     MagnetArray(3)


    %-----------------------------------------------------------------------------------------------
    % Directory/File preparation and creation
    %----------------------------------------
    run.save.updateDataSaves();
    
    if ~run.save.FSAVE; return; end
    
    iteration = run.time.iteration;
    saveState = uint8(0);
    saveState = bitset(saveState, 1, iteration == 0); % Beginning of a run
    saveState = bitset(saveState, 2, run.save.done); % End of run

    %-----------------------------------------------------------------------------------------------
    % Save Data to File
    %------------------
    if run.save.saveData
        
        %--- Preparation for save using lab 1 ---%
        if labindex == 1
            switch saveState
                case 2;            
                    fileSuffix         = 'FINAL'; 
                    
                case 1;            
                    fileSuffix         = 'START';
                    run.save.firstSave = false; % Update for later saves.
                    
                otherwise;        
                    fileSuffix         = run.paths.iterationToString(iteration);
            end
        
            sl.gamma  = run.GAMMA;
            sl.time   = run.time.toStruct();
            sl.about  = run.about;
            sl.ver    = run.version;
            sl.iter   = iteration;
            
            if run.treadmill.ACTIVE,    sl.tread = run.treadmill.toStruct(); end
            
        end

        for i=1:size(run.save.SLICE,1) % For each slice type
            switch (i)
                case {1 2 3};    % 1D Slices
                    if ( ~run.save.save1DData ); continue;        else sliceDim = '1D';     end
                case {4 5 6};    % 2D Slices
                    if ( ~run.save.save2DData ); continue;        else sliceDim = '2D';     end
                case 7;            % 3D Slice
                    if ( ~run.save.save3DData ); continue;        else sliceDim = '3D';     end
                case 8;            % Custom Slice
                    if ( ~run.save.saveCustomData ); continue;    else sliceDim = 'Cust';   end
            end

            %--- Save slice if active ---%
            if (run.save.ACTIVE(i)) % Check if active
                sl.mass = run.save.getSaveSlice(mass.array, i);
                sl.momX = run.save.getSaveSlice(mom(1).array, i);
                sl.momY = run.save.getSaveSlice(mom(2).array, i);
                sl.momZ = run.save.getSaveSlice(mom(3).array, i);
                sl.ener = run.save.getSaveSlice(ener.array, i);
                
                if (run.magnet.ACTIVE && ~isempty(mag(1).array))
                    sl.magX = run.save.getSaveSlice(mag(1).array, i);
                    sl.magY = run.save.getSaveSlice(mag(2).array, i);
                    sl.magZ = run.save.getSaveSlice(mag(3).array, i);
                else sl.magZ = []; sl.magY = []; sl.magX = [];
                end
                
                if (run.gravity.ACTIVE && numel(grav.array) > 1)
                    sl.grav = run.save.getSaveSlice(grav.array, i);
                else sl.grav = []; 
                end
                
                %--- Slice DGRID to match arrays if necessary ---%
                sl.dGrid = cell(1,3);
                for n=1:3
                    if numel(run.DGRID{n}) > 1
                        sl.dGrid{n} = run.save.getSaveSlice(run.DGRID{n}, i);
                    else
                        sl.dGrid{n} = run.DGRID{n};
                    end
                end
                
                %--- Save only using lab 1 ---%
                if labindex == 1
                    sl.dim = sliceDim;

                    fileName = [run.paths.save, '/', sliceDim, '_', run.save.SLICELABELS{i}, ...
                                '_', fileSuffix];
                    sliceName = strcat('sx_', run.save.SLICELABELS{i}, '_', fileSuffix);
                    if ~isvarname(sliceName); sliceName = genvarname(sliceName); end
                    
                    try
                        eval( [sliceName ' = sl;'] );
                        eval( ['save(''' fileName ''',''' sliceName ''');' ] );
                    catch MERR %#ok<NASGU>
                        fprintf('Unable to save. Skipping');
                    end
                end
                labBarrier(); % Block all labs until lab 1 is ready to proceed.
            end
        end
    end

    %-----------------------------------------------------------------------------------------------
    % Save Info File
    %---------------
    if (saveState && labindex == 1)
                
        %-------------------------------------------------------------------------------------------
        % Store Array Information
        %------------------------
        if run.save.done
            fid = fopen([run.paths.save '/arrays-end.log'],'a');
            fprintf(fid,'---++ Arrays at Complete:');
        else
            fid = fopen([run.paths.save '/arrays-ini.log'],'a');
            fprintf(fid,'---++ Arrays at Creation:');
        end

        fprintf(fid, '\n\n---+++ Run Manager\n%s', ImogenRecord.valueToString(run, {'parent'}) );
        fprintf(fid, '\n\n---+++ Mass\n%s',        ImogenRecord.valueToString(mass) );
        fprintf(fid, '\n\n---+++ Momentum\n%s',    ImogenRecord.valueToString(mom) );
        fprintf(fid, '\n\n---+++ Energy\n%s',      ImogenRecord.valueToString(ener) );
        fprintf(fid, '\n\n---+++ Magnet\n%s',      ImogenRecord.valueToString(mag) );
        fprintf(fid, '\n\n---+++ Gravity\n%s',     ImogenRecord.valueToString(grav) );
        fclose(fid);
        
            
        %---------------------------------------------------------------------------------------
        % Create runInfo file
        %--------------------

        %--- Get computer and account information ---%
        [~, host] = system('hostname'); host = deblank(host);
        [~, user] = system('whoami');   user = deblank(user);

        %--- Determine the amount of time elapsed during the run ---%
        endVec = clock; %Get current time
        if run.save.done
            dateNumStart    = datenum(run.time.startTime); dateNumEnd = datenum(endVec);
            dDateNum        = dateNumEnd - dateNumStart;
            dDays           = floor(dDateNum);
            dHours          = floor( 24 * (dDateNum - dDays) );
            dMinutes        = 60 * (24 * (dDateNum - dDays) - dHours);                      
            delTimeStr      = sprintf('[%g days][%g hours][%g minutes]',dDays, dHours, dMinutes);    
            logName         = '/runInfo.log';
        else logName = '/start.log';
        end
        
        fid = fopen([run.paths.save, logName],'a');

        fprintf(fid, sprintf( '%s\n', run.about ));
        fprintf(fid,'\n---++ Run Details\n');
        fprintf(fid,'   * *Save code:* %s\n', run.paths.saveFolder);
        fprintf(fid,'   * *Imogen code:* v.%s\n', num2str(run.version, '%0.2g'));
        fprintf(fid,'   * Started at: %s\n', datestr(run.time.startTime,'HH:MM on mmm dd, ''yy'));
        if run.save.done
            fprintf(fid,'   * Ended at:   %s\n', datestr(endVec  ,'HH:MM on mmm dd, ''yy'));
            fprintf(fid,'   * Elapsed Run time: %s\n', delTimeStr);
        end
        
        fprintf(fid,'   * Run on: %s\n   * by user: %s\n', host, user);
        fprintf(fid,'   * Using Matlab %s (v. %s)\n', run.matlab.Release, run.matlab.Version);

        fprintf(fid, '\n---++ Run information');
        fprintf(fid, sprintf([run.notes '\n']));
        
        fprintf(fid, '\n---++ Warnings');
        fprintf(fid, sprintf([run.warnings '\n']));
        
        fprintf(fid,['\n<!-- Initialization Parameters -->\n%%TWISTY{prefix="<h2>" mode="div" ' ...
                     'link="Initialization" showimgleft=', ...
                     '"%%ICONURLPATH{toggleopen}%%" hideimgleft="%%ICONURLPATH{toggleclose}%%" ', ...
                     'class="twikiHelp" suffix="</h2>"}%%\n']);
        fprintf(fid, sprintf(run.iniInfo));
        fprintf(fid, '\n%%ENDTWISTY%%\n');
        
        fprintf(fid,['\n<!-- Settings Parameters -->\n%%TWISTY{prefix="<h2>" mode="div" ' ...
                     'link="Settings" showimgleft=', ...
                     '"%%ICONURLPATH{toggleopen}%%" hideimgleft="%%ICONURLPATH{toggleclose}%%" ', ...
                     'class="twikiHelp" suffix="</h2>"}%%\n']);
        if run.save.done
            fprintf(fid, sprintf( '   * Iterations: %d\n'      , iteration) );
            fprintf(fid, sprintf( '   * Simulated Time: %g\n'  , run.time.time) );
        end

        %--- Print entire infoArray ---%
        riLen = size(run.info);
        for i=1:riLen(1)
            if (~isempty(run.info{i,1}))
                if ~ischar(run.info{i,2}), run.info{i,2} = mat2str(run.info{i,2}); end
                fprintf(fid, sprintf( '   * %s: %s\n',run.info{i,1},run.info{i,2}) );
            end
        end

        %--- Print Gravity solver information ---%
        if ~isempty(run.gravity.info), fprintf(fid, ['\n---++ Gravity\n' run.gravity.info]);
        else fprintf(fid, '   * Gravity solver performed without error.');
        end
        fprintf(fid, '\n%%ENDTWISTY%%\n');
        fclose(fid);

        %--- Notify of code of finish ---%
        if run.save.done
            fprintf('Run completed at %s.\n',datestr(endVec ,'HH:MM on mmm-dd-yy'));
            disp(strcat('Results files saved to directory: ', run.paths.save));
        end
    end
    labBarrier(); % Block all labs until lab 1 has finished saving

    % Save images
    run.image.imageSaveHandler(mass, mom, ener, mag, grav);
    
end
