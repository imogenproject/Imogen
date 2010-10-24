function makeVectorFile(filebase, stepno, varwidth, Mx, My, Mz, vardesc)

fname = sprintf(['%s.%0' num2str(varwidth, '%i') 'i'], filebase, stepno);
VEC = fopen(fname, 'w');

% Write description & 'part'
charstr = char(32*ones([160 1]));
charstr(1:length(vardesc)) = vardesc;
charstr(81:84)             = 'part';
fwrite(VEC, charstr, 'char*1');

% Write part number
fwrite(VEC, [1], 'int');

% write 'block'
charstr = char(32*ones([80 1]));
charstr(1:5) = 'block';
fwrite(VEC, charstr, 'char*1');

% write all x-direction, ydirection and zdirection vectors
fwrite(VEC, reshape(Mx, [prod(size(Mx)) 1]), 'float');
fwrite(VEC, reshape(My, [prod(size(My)) 1]), 'float');
fwrite(VEC, reshape(Mz, [prod(size(Mz)) 1]), 'float');

fclose(VEC);

end
