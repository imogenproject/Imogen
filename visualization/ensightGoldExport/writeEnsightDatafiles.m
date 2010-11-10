function exportEnsightDatafiles(basename, frameNo, frame)
% basename: output filename base
% frameNo:  frame number
% frame: Imogen sx_... savefile structure

makeEnsightScalarFile(sprintf('%s.mass.%04i', basename, frameNo), frame.mass, 'mass');
makeEnsightScalarFile(sprintf('%s.ener.%04i', basename, frameNo), frame.ener, 'energy');

if ~isempty(frame.grav)
    makeEnsightScalarFile(sprintf('%s.grav.%04i', basename, frameNo), frame.grav, 'gravity');
end

makeEnsightVectorFile(sprintf('%s.mom.%04i', basename, frameNo), ...
                      frame.momX, frame.momY, frame.momZ, 'momentum');

if ~isempty(frame.magX)
    makeEnsightVectorFile(sprintf('%s.mag.%04i', basename, frameNo), ...
                          frame.magX, frame.magY, frame.magZ, 'magnetic_field');
end


end
