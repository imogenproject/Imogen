function magCell = getMagnetCell_Util(mag)


	magCell = zeros(size(mag));
	
	switch ndims(magCell)
		case 2
			magCell = 0.5 * ( mag + circ_shift(mag,2,1) );
			
		case 3
			for i=1:2
				magCell(i,:,:) = 0.5 * ( squeeze(mag(i,:,:)) ...
											+ circ_shift( squeeze(mag(i,:,:)),i,1) );
			end
			
		case 4
			for i=1:3
				magCell(i,:,:,:) = 0.5 * ( squeeze(mag(i,:,:,:)) ...
											+ circ_shift( squeeze(mag(i,:,:,:)),i,1) );
			end
    
end