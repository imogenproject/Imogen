function [resArray,sampledArray] = analyze_SedovTaylor(array,center,bounds,cutoffDist,sampleRate)
% Analyzes array results for the Sedov-Taylor blast wave by collecting the spherically symmetric 
% profiles and displaying them in an overlayed fashion to verify the Cartesian grid error for a
% spherical representation.
%
%>> array			array to analyze											double	[nx ny nz]
%>> center			grid point that represents the center of the blast wave		double	[x y z]
%>> bounds			
%>> cutoffDist
%>> sampleRate
%<< resArray		resulting array from the sampled space
%<< sampledArray	


    if (nargin < 1 || isempty(array)), error('Imogen:DataInputError','No array specified. Operation aborted.'); end
    nDim = ndims(array);
    grid = size(array);
    resArray = zeros([5, prod(grid)]);
    
    if (nargin < 5 || isempty(sampleRate)), sampleRate = 0.01; end
    if (nargin < 4 || isempty(cutoffDist)), cutoffDist = 1.0; end
    if (nargin < 3 || isempty(bounds)), bounds = [minFinderND(array), maxFinderND(array)]; end
    if (nargin < 2 || isempty(center)), center = round(grid / 2.0); end

    maxDist = sqrt(sum((grid - center).^2,2));
    
    if (nDim == 3)
        n = 1;
        for i=1:grid(1)
            for j=1:grid(2)
                for k=1:grid(3)
                    if ( (array(i,j,k) > bounds(1)) && (array(i,j,k) < bounds(2)) )
                        index = [i j k] - center;
                        distVal = sqrt(sum(index.*index,2))/maxDist;
                        if (distVal < cutoffDist)
                            resArray(1:2,n) = [distVal array(i,j,k)];
                            resArray(3:5,n) = abs(index) / sqrt(sum(index.*index,2));
                            n = n + 1;
                        end
                    end
                end
            end
        end
    end
    resArray = resArray(:,1:(n-1));
    
    sampleSize = floor(min(sampleRate * prod(grid),0.75*(n-1)));
    rVals{1,2} = ceil((n-1) * rand(1,sampleSize));
    rVals{1,1} = 1:5;
    sampledArray = resArray(rVals{:});
    
end