function makeVectorFile(filebase, Mx, My, Mz, vardesc)

VEC = fopen(filebase, 'w');

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
fwrite(VEC, reshape(single(Mx), [prod(size(Mx)) 1]), 'float');
fwrite(VEC, reshape(single(My), [prod(size(My)) 1]), 'float');
fwrite(VEC, reshape(single(Mz), [prod(size(Mz)) 1]), 'float');

fclose(VEC);

end
