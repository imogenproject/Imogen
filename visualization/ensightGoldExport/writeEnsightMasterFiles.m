function writeEnsightMasterFiles(basename, range, frame, timeNormalization)

CASE = fopen([basename '.case'], 'w');

% prepare format info
fprintf(CASE, 'FORMAT\ntype: ensight gold\n');

% prepare geometry info
fprintf(CASE, '\nGEOMETRY\n');
fprintf(CASE, 'model: 1 %s.geom\n', basename);
makeEnsightGeometryFile(frame, basename);

% prepare variables
fprintf(CASE, '\nVARIABLE\n');

fprintf(CASE, 'scalar per node: 1 mass %s.mass.****\n', basename);
fprintf(CASE, 'scalar per node: 1 energy %s.ener.****\n', basename);
if ~isempty(frame.grav)
    fprintf(CASE, 'scalar per node: 1 grav_potential %s.grav.****\n', basename);
end

fprintf(CASE, 'vector per node: 1 momentum %s.mom.****\n', basename);

if ~isempty(frame.magX)
    fprintf(CASE, 'vector per node: 1 magnet %s.mag.****\n', basename);
end

fprintf(CASE, '\nTIME\n');
fprintf(CASE, 'time set:              1 time_data\n');
fprintf(CASE, 'number of steps:       %i\n', numel(range));
fprintf(CASE, 'filename start number: 0\n');
fprintf(CASE, 'filename increment:    1\n');

nwritten = fprintf(CASE, 'time values: ');
for q = 1:numel(range)
    nwritten = nwritten + fprintf(CASE,'%5.5g ', sum(frame.time.history(1:range(q)))/timeNormalization);
    if nwritten > 72; fprintf(CASE, '\n'); nwritten = 0; end
end
fprintf(CASE, '\n');
fclose(CASE);

% Don't want this unless doing parallel data apparently.
%fprintf(CASE, '\n\nFILE\n');
%fprintf(CASE, 'file set: 1\n');
%fprintf(CASE, 'number of steps: %i\n', numel(range)); 

end
