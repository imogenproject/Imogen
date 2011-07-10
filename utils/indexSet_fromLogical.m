function result = indexSet(logic)
% Given the dimensions in dims and sets of x, y and z coordinates,
% This function takes the volume spanned by the cartesian of xset x yset x zset
% and returns [linear x y z] Nx4 matrix corresponding to the linear and x y z
% indicies corresponding to the selected points

dims = size(logic);
if numel(dims) == 2; dims(3) = 1; end

[u v w] = ndgrid(1:size(logic,1), 1:size(logic,2), 1:size(logic,3) );

%u = u(:);
%v = v(:);
%w = w(:);

result = [];

j = 1;

for i = 1:numel(logic)
    if(logic(i) == 1) 
        lind = 1+(u(i)-1)+(v(i)-1)*dims(1) + (w(i)-1)*dims(1)*dims(2);
        result(j,:) = [lind u(i) v(i) w(i)];
        j=j+1;
    end
end

end 
