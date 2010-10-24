function [lowerFactorization upperFactorization] = poissonBlockILU(M, dropTolerance, blockSize, matSize)
% This function computes a block incomplete LU (Block ILU) factorization of the coefficient matrix
% describing the Poisson equation. This matrix is real, symmetric and regular
% Matlab appears to be very inefficient at inserting large blocks of data into already-existing
% sparse matrices so we will store large vectors of I, J, and K values (matlab: 'doc sparse')
% and then create the sparse matrix once.
%
%>> M                   Large sparse matrix describing discretized Poisson equation     sparse
%>> dropTolerance       Factorization drop tolerance. For explanation of this           double
%                           see property see help for the luinc command 'help luinc',
%                           which uses the property.
%>> blockSize           Size of the blocks; This must be chosen to be an integer        int
%                           fraction of the size of M or bad, BAD things will 
%                           happen.
%<< lowerFactorization  Complete lower factorization matrix result                      sparse
%<< upperFactorization  Complete upper factorization matrix result                      sparse

    %--- Compute the LU factorization of one block ---%
    %       M must be a regular matrix (of the type associated with a discrete partial difference
    %       operator like the discrete Laplacian for example) and blockSize must match the size of 
    %       a repeating feature This will make every block identical and let us get away with only 
    %       one LU factorization.
    fprintf('   Computing LU factor... ');
    %matSize     = [prod(matsize) prod(matsize)];
    imax        = matSize(1)/blockSize;
    blockIdx    = 0;
    blksub      = M(blockSize*blockIdx + 1:blockSize*blockIdx + blockSize , ...
                    blockSize*blockIdx + 1:blockSize*blockIdx + blockSize);
    [lf uf]     = luinc(blksub, dropTolerance);

    %--- Begin replicating the resulting factor down the diagonal ---%
    %        6 vectors hold components for 2 sparse matrices
    %        Define here to be persistent through the loop below
    %        Then extract the nonzero indices and entries of the l/u factors
    bigIL = []; bigIU = [];
    bigJL = []; bigJU = [];
    bigKL = []; bigKU = [];

    %--- Pick out the nonzero parts of the block ---%
    nonzerolower    = (lf ~= 0);
    nonzeroupper    = (uf ~= 0);

    [IL JL]         = find(nonzerolower);
    KL              = lf(nonzerolower);

    [IU JU]         = find(nonzeroupper);
    KU              = uf(nonzeroupper);

    %--- Iteratively append factor values onto big vectors and move them down the diagonal ---%
    fprintf('Writing matrix (%i *): ', imax);
    for i = 0:(imax-1)
        bigIL = [bigIL IL]; bigJL = [bigJL JL]; bigKL = [bigKL KL];
        bigIU = [bigIU IU]; bigJU = [bigJU JU]; bigKU = [bigKU KU];

        IL = IL + blockSize; JL = JL + blockSize;
        IU = IU + blockSize; JU = JU + blockSize;
        fprintf('*');
        if mod(i, 50) == 49; fprintf('\n'); end;
    end

    lowerFactorization = sparse(bigIL, bigJL, bigKL, matSize(1), matSize(2));
    upperFactorization = sparse(bigIU, bigJU, bigKU, matSize(1), matSize(2));
    fprintf('\n');

end
