#include "mex.h"
#include "matrix.h"

#include "math.h"
#include "time.h"

int lvlMax;
double rconst;

double sumLevel(double **rhos, double **poss, double *r, int lvl, int *lvlSize, int *bounds);

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
/* Requires 5 arguments in this order:
 1. Cell array of hierarhcially coarsened density
 2. Cell array of hierarchially coarsened mass centers
 3. Realspace locations to sum at (9x1 vector)
 4. Coarsening algorithm's constant
%% NOT ANYMORE 5. Mean cell distance from center (1/8 length of line to opposite corners of cell)
*/

  /* Check for proper number of arguments. */
if(nrhs!=4) {
	mexErrMsgTxt("four inputs required.");
} else if(nlhs>1) {
	mexErrMsgTxt("One return value required.");
}

/* Get number of argument and allocate memory */
int nLevels = mxGetNumberOfElements(prhs[0]);

int *arraydims = malloc(3*nLevels * sizeof(int));
double **rhoptr = malloc(nLevels * sizeof(double*));
double **posptr = malloc(nLevels * sizeof(double*));

rconst = *mxGetPr(prhs[3]);

int j;
/* Build arrays of the size of every level and *s to them */
for(j = 0; j < nLevels; j++) {
	mxArray *rhos = mxGetCell(prhs[0], j);
	mxArray *poss = mxGetCell(prhs[1], j);

	mwSize d = mxGetNumberOfDimensions(rhos);
	mwSize *s = mxGetDimensions(rhos);

	if(d == 2) {
		/* 2-D */
		
		arraydims[3*j] = (int)s[0];
		arraydims[3*j+1]=(int)s[1];
		arraydims[3*j+2]= 1;
		} else {
		/* 3-D, it better not be higher-dimensional */
		arraydims[3*j] = (int)s[0];
		arraydims[3*j+1]=(int)s[1];
		arraydims[3*j+2]=(int)s[2];
		}

	rhoptr[j] = mxGetPr(rhos);
	posptr[j] = mxGetPr(poss);
	}

/*double *cmr = mxGetPr(prhs[4]); */

/* At some point this will get broken up into a parallel operation, then a CUDA operation. */
int sumvol[6] = {0, 0, 0, arraydims[0], arraydims[1], arraydims[2]};

/* 9x1 - [ [startxyz] [endxyz] [npointsxyz] ] */
double *vin = mxGetPr(prhs[2]);

double r0[3] = {vin[0], vin[1], vin[2]};
double dr[3] = { (vin[3]-vin[0])/vin[6], (vin[4]-vin[1])/vin[7], (vin[5]-vin[2])/vin[8] };

for(j = 0; j < 3; j++) if(vin[6+j] == 0) dr[j] = 0; /* Prevent divide by zero */

int ctx, cty, ctz;

int maxid[3] = { (int)vin[6]+1, (int)vin[7]+1, (int)vin[8]+1 };

mwSize retdims[3] = { maxid[0], maxid[1], maxid[2] };
plhs[0] = mxCreateNumericArray(3, retdims, mxDOUBLE_CLASS, mxREAL);

double *phi = mxGetPr(plhs[0]);
double rp[3];

lvlMax = nLevels;

for(ctx = 0; ctx < maxid[0]; ctx++)
for(cty = 0; cty < maxid[1]; cty++)
for(ctz = 0; ctz < maxid[2]; ctz++) {
	rp[0] = r0[0] + ctx*dr[0];
	rp[1] = r0[1] + cty*dr[1];
	rp[2] = r0[2] + ctz*dr[2];

	phi[ctx + maxid[0]*cty + maxid[0]*maxid[1]*ctz] = sumLevel(rhoptr, posptr, rp, 0, arraydims, sumvol);
}

free(arraydims); free(rhoptr); free(posptr);

return;
}

double metric(double *a, double *b)
{
/* This is correct for the three-dimensional 1/r potential.
 Can be modified to return 1/ln(r) for a line potential or 1/r for a plane */
/*printf("	Pos: <%3f %3f %3f>\n", b[0], b[1], b[2]);*/

return sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]));
}

double sumLevel(double **rhos, double **poss, double *r, int lvl, int *lvlSize, int *bounds)
{
double dphi = 0;

int i, j, k; /* Array indices */
int *size = &lvlSize[3*lvl];

double *rho = rhos[lvl];
double *pos = poss[lvl];

int poslidx = 0;
int rholidx = 0;
double rad;

int rid[3]; 


for(i = bounds[0]; i < bounds[3]; i++)
for(j = bounds[1]; j < bounds[4]; j++)
for(k = bounds[2]; k < bounds[5]; k++) {
	/* Examine every cell at this level */
 	rholidx = (i + size[0]*j + (size[0]*size[1])*k);
	poslidx = 3*rholidx;

	rad = metric(r, &pos[poslidx]);

/*printf("%i %i %i %i: %g ", i+1, j+1, k+1, lvl+1, rad);*/

	if(rad > (rconst * (double)(1 << (lvlMax - lvl)))) {
		/* If we're far enough away, sum this cell */
		dphi -= rho[rholidx] / rad;
/*printf("*\n");*/
		} else {
/*printf("\n");*/
		/* Too close, recurse if possible */
		if (lvl < (lvlMax-1)) {
			int newbounds[6] = { 2*i, 2*j, 2*k, 2*i+2, 2*j+2, 2*k+2 };
	
			if(newbounds[3] > lvlSize[3*lvl+3]) newbounds[3]--;
			if(newbounds[4] > lvlSize[3*lvl+4]) newbounds[4]--;
			if(newbounds[5] > lvlSize[3*lvl+5]) newbounds[5]--;

			dphi += sumLevel(rhos, poss, r, lvl+1, lvlSize, newbounds);
			} else {
			/* If we're not in the same cell we're summing, use finest level available.
			 If we are, use the mean distance of the cell's contents from the center */
			if(rad > 1e-6) { dphi -= rho[rholidx]/rad; } /*else { dphi -= rho[rholidx] / cmr; } */
			}
		}
	}

return dphi;

}


