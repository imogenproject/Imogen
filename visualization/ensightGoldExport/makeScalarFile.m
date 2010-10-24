function makeScalarFile(filebase, stepno, varwidth, M, vardesc)

fname = sprintf(['%s.%0' num2str(varwidth, '%i') 'i'], filebase, stepno);
SCAL = fopen(fname, 'w');

% Write description & 'part'
charstr = char(32*ones([160 1]));
charstr(1:length(vardesc)) = vardesc;
charstr(81:84)             = 'part';
fwrite(SCAL, charstr, 'char*1');

% Write part number
fwrite(SCAL, [1], 'int');

% write 'block'
charstr = char(32*ones([80 1]));
charstr(1:5) = 'block';
fwrite(SCAL, charstr, 'char*1');

% write all scalars in array M
fwrite(SCAL, reshape(M, [prod(size(M)) 1]), 'float');

fclose(SCAL);

end
