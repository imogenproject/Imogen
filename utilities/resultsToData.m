function data = resultsToData(result)
% function to convert Imogen results into a structure containing data stored as classes. This
% is useful if you want to use functions and operators that require object input instead of arrays.

    data.run  = fauxImogenManager(result.dGrid, size(result.mass));
    data.mass = FluidArray(ENUM.SCALAR, ENUM.MASS, result.mass, data.run, []);
    data.ener = FluidArray(ENUM.SCALAR, ENUM.ENER, result.ener, data.run, []);
    
    data.grav = GravityArray(ENUM.GRAV, data.run, []);
    if isempty(result.grav);        data.grav.array = zeros(size(result.mass));
    else                            data.grav.array = data.grav;
    end

    data.mom  = FluidArray.empty(3,0);
    data.mag  = MagnetArray.empty(3,0);

    fields = {'X','Y','Z'};
    for i=1:3
        data.mom(i) = FluidArray(ENUM.VECTOR(i), ENUM.MOM, ...
                                 result.(['mom',fields{i}]), data.run, []);
                             
         if isempty(result.(['mag' fields{i}]))
             data.mag(i) = MagnetArray(ENUM.VECTOR(i), ENUM.MAG, ...
                                       zeros(size(result.mass)), data.run, []);
         else
             data.mag(i) = MagnetArray(ENUM.VECTOR(i), ENUM.MAG, ...
                                       result.(['mag',fields{i}]), data.run, []);
         end                    
    end
    
end
