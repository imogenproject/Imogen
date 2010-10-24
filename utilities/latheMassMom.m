function [mass momX momY] = latheMassMom(rhoSlice, momSlice, gridres, h, isAngMom)
% Assumes zhat is the cylindrical axis

mass = zeros(gridres);
momX = zeros(gridres);
momY = zeros(gridres);

radVector = (1:size(rhoSlice,1))'*h;

for zct = 1:gridres(3)
    if isAngMom
        momSlice(:,zct) = momSlice(:,zct) ./ radVector;
        momSlice(isnan(momSlice)) = 0; % Fix 1/0 problems if they occur
    end
    [mass(:,:,zct) momX(:,:,zct) momY(:,:,zct)] = cyl2rect(radVector, rhoSlice(:,zct), momSlice(:,zct), gridres(1)/2, h);
end

end
