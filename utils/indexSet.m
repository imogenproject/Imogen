function result = indexSet(dims, xset, yset, zset)
% Given the dimensions in dims and sets of x, y and z coordinates,
% This function takes the volume spanned by the cartesian of xset x yset x zset
% and returns [linear x y z] Nx4 matrix corresponding to the linear and x y z
% indicies corresponding to the selected points

[u v w] = ndgrid(xset, yset, zset);

u = u(:);
v = v(:);
w = w(:);

lind = 1+(u-1)+(v-1)*dims(1) + (w-1)*dims(1)*dims(2);

result = [lind u v w];

end 
