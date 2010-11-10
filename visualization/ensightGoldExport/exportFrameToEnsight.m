function exportToEnsightTest(frame)
% FIXME: Add stuff to support not loading everything first

basename = input('Basename: ','s');

exportEnsightDatafiles(basename, 0, frame);
writeEnsightMasterFiles(basename, numel(frame.time.history), frame, 1);

end
