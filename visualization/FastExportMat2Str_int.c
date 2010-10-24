#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
 int i, j, m, n;
 
 
 double *data1, *data2;
 
 if (nrhs != nlhs)
 mexErrMsgTxt("The number of input and output arguments must be the same.");
 
 char *mainstr = NULL;
 
 for (i = 0; i < nrhs; i++) 
   {
    /* Find the dimensions of the data */
    m = mxGetM(prhs[i]);
    n = mxGetN(prhs[i]);
 
    int N = m*n*7; /* Max amount of storage needed for our string */
    
    mainstr = realloc(mainstr, N); /* Allocate storage */

/* printf("%i: %x\n", N, mainstr); fflush(stdout); */

    mainstr[0] = '[';
    
    data1 = mxGetPr(prhs[i]); /* Get pointer to given input string */
    int ppos = 1;
    for (j = 0; j < m*n; j++)
    {
    ppos += sprintf(&mainstr[ppos],  "%i ", (int)(data1[j]) );

    /* Grow if we get too large */
    if( (N - ppos) < 100) { N = N + N/4; mainstr = realloc(mainstr, N); }
    }
    sprintf(&mainstr[ppos], "]");
    
/*printf("%i\n", ppos); fflush(stdout);*/

    /* Create an mxArray for the output data */
    plhs[i] = mxCreateString(mainstr);

/*printf("Created mxString\n"); fflush(stdout);*/
    
   }

free(mainstr);
 
}
