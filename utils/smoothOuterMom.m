function s = smoothOuterMom(r, dist)

s = 1.0*r;
size(r)
for x = 1:size(r, 1)
for y = 1:size(r, 2)
  z = size(r,3);
  tripped = 0;

  while (tripped < dist) && (z > 1)
    z = z - 1;
    if abs(r(x,y,z)) > 0
      s(x,y,z) = r(x,y,z)*(tripped+.5) / dist;
      tripped = tripped + 1;
    end
  end

  z = 1;
  tripped = 0;
  while (tripped < dist) && (z < (size(r,3)-1 ))
    z = z + 1;
    if abs(r(x,y,z)) > 0
      s(x,y,z) = r(x,y,z)*(tripped+.5) / dist;
      tripped = tripped + 1;
    end
  end

  
end
end

end 
