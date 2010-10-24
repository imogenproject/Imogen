function imogenParallelLoad(runFile, logFile, threads)
% This script is the command and control manager for command line run Imogen scripts in a parallel 
% environment. It is designed to simplify the run script syntax and allow for greater extensibility.
% 
%>> runFile		run file function name to execute								str
%>> logFile		log file name for writing output information					str
%>> threads		number of parallel threads for distributed computing			int

	pjob = createParallelJob();
	set(pjob,'MinimumNumberOfWorkers',threads);
	set(pjob,'MaximumNumberOfWorkers',threads);
	createTask(pjob, 'imogenLoad', 0, {runFile, logFile});
	assignin('base','pjob',pjob);
	submit(pjob);
	waitForState(pjob,'finished');
	fprintf('\n------ Job information ------\n');
	pjob.display();
	pjob.Tasks.display();
	for i=1:length(pjob.Tasks)
		fprintf('\n------- Task %g information ------\n',i);
		pjob.Tasks(i).display;
	end
	destroy(pjob);
end