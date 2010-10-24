function exportLegacyVTK(frame, filename)

fid = fopen([filename '.vts'],'w','b');

Nx = size(frame.mass, 1);
Ny = size(frame.mass, 2);
Nz = size(frame.mass, 3);
Ntot = Nx * Ny * Nz;

pvsc = 0;
%fwrite(fid, string)
fwrite(fid, sprintf('# vtk DataFile Version 2.0\nTesting - imogen frame save export\nBINARY\n'));

if (numel(frame.dGrid{1}) > 1) || (numel(frame.dGrid{2} > 1)) || (numel(frame.dGrid{3}) > 1)
    % Print dx, dy, dz vectors
    fwrite(fid, sprintf('DATASET RECTILINEAR_GRID\nDIMENSIONS %i %i %i\nX_COORDINATES %i float\n', Nx+pvsc, Ny+pvsc, Nz+pvsc, Nx+pvsc));
    
    if numel(frame.dGrid{1}) > 1
        q = cumsum(frame.dGrid{1}(:,1,1)); q = q - q(1);
    else
        q = (0:(Nx-1))*frame.dGrid{1};
    end
    fwrite(fid, q, 'float');
    
    fwrite(fid, sprintf('\nY_COORDINATES %i float\n', Ny+pvsc));
    if numel(frame.dGrid{2}) > 1
       q = cumsum(squeeze(frame.dGrid{2}(1,:,1))); q = q - q(1);
    else
        q = (0:(Ny-1))*frame.dGrid{2};
    end
    fwrite(fid, q, 'float');
    
    fwrite(fid, sprintf('\nZ_COORDINATES %i float\n', Nz+pvsc));
    if numel(frame.dGrid{3}) > 1
        q = cumsum(squeeze(frame.dGrid{3}(1,1,:))); q = q - q(1);
    else
        q = (0:(Nz-1))*frame.dGrid{3};
    end
    fwrite(fid, q, 'float');
    fwrite(fid, sprintf('\n'));
else
    % Uniform rectilinear coordinates
    fwrite(fid, sprintf('DATASET STRUCTURED_POINTS\nDIMENSIONS %i %i %i\nORIGIN 0 0 0\nSPACING %g %g %g\n', Nx+pvsc, Ny+pvsc, Nz+pvsc, frame.dGrid{1}, frame.dGrid{2}, frame.dGrid{3}));
end

fprintf('\n\tWriting mass... ');
% Mass density
%dataStr = FastExportMat2Str( reshape(frame.mass, [1 Ntot] ) );
fwrite(fid, sprintf('POINT_DATA %i\nSCALARS mass float 1\nLOOKUP_TABLE default\n', Ntot),'char*1');
fwrite(fid, reshape(frame.mass, [1 Ntot]), 'float');

fprintf('Done.\n\tWriting energy... ');
% Energy density
%datastr = FastExportMat2Str( reshape(frame.ener, [1 Ntot] ) );
fwrite(fid, sprintf('SCALARS energy float 1\nLOOKUP_TABLE default\n'));
fwrite(fid, reshape(frame.ener, [1 Ntot]), 'float');

fprintf('Done.\n\tWriting gravitational potential... ');

if numel(frame.grav) == numel(frame.mass)
%	datastr = FastExportMat2Str( reshape(frame.grav, [1 Ntot]) );
	fwrite(fid, sprintf('SCALARS gravpot float 1\nLOOKUP_TABLE default\n'));
	fwrite(fid, reshape(frame.grav, [1 Ntot]), 'float');

	fprintf('Done.\n\tWriting momentum... ');
else
	fprintf('Empty.\n\tWriting momentum... ');
end

% Momentum
momvec(1:3:(3*Ntot)) = reshape(frame.momX, [1 Ntot]);
momvec(2:3:(3*Ntot)) = reshape(frame.momY, [1 Ntot]);
momvec(3:3:(3*Ntot)) = reshape(frame.momZ, [1 Ntot]);

%dataStr = FastExportMat2Str( reshape(momvec, [1 3*Ntot]) );
fwrite(fid, sprintf('VECTORS momentum float\n'));
fwrite(fid, momvec, 'float');

clear momvec;
fprintf('Done.\n\tWriting magnetic field... ');

if numel(frame.magX) == numel(frame.mass)
	magvec(1:3:(3*Ntot)) = reshape(frame.magX, [1 Ntot]);
	magvec(2:3:(3*Ntot)) = reshape(frame.magY, [1 Ntot]);
	magvec(3:3:(3*Ntot)) = reshape(frame.magZ, [1 Ntot]);
	
%	dataStr = FastExportMat2Str( reshape(magvec, [1 3*Ntot]) );
	fwrite(fid, sprintf('VECTORS magnetic float\n'));
	fwrite(fid, magvec, 'float');
	
	fprintf('Done.\n');
	clear magvec;
else
	fprintf('Empty.\n');
end



end
