function exportToEnsightTest(frame)
% FIXME: Add stuff to support not loading everything first

basename = input('Basename: ','s');
CASE = fopen([basename '.case'], 'w');

% prepare format info
fprintf(CASE, 'FORMAT\ntype: ensight gold\n');

% prepare geometry info
fprintf(CASE, '\nGEOMETRY\n');
fprintf(CASE, 'model:		%s.geom\n', basename);
makeGeometryFile(frame, basename);

% prepare variables
fprintf(CASE, '\nVARIABLE\n');

fprintf(CASE, 'scalar per node:		mass %s.mass.**\n', basename);
makeScalarFile([basename '.mass'], 0, 2, frame.mass, 'mass');

fprintf(CASE, 'scalar per node: 	energy %s.ener.**\n', basename);
makeScalarFile([basename '.ener'], 0, 2, frame.ener, 'energy');

if ~isempty(frame.grav)
    fprintf(CASE, 'scalar per node:             grav_potential %s.grav.**\n', basename);
    makeScalarFile([basename '.grav'], 0, 2, frame.grav, 'gravity');
end

fprintf(CASE, 'vector per node:		momentum %s.mom.**\n', basename);
makeVectorFile([basename '.mom'], 0, 2, frame.momX, frame.momY, frame.momZ, 'momentum');

if ~isempty(frame.magX)
    fprintf(CASE, 'vector per node:		magnet %s.mag.**\n', basename);
    makeVectorFile([basename '.mag'], 0, 2, frame.magX, frame.magY, frame.magZ, 'magnetic_field');
end

% FIXME: We need time information and the TIME field here

end
