#include "mex.h"
#include <cmath>
void mexFunction(int nlhs, mxArray *plhs[],	int nrhs, const mxArray *prhs[])
{
	double *X;
	double *Y;
	bool *shield;
	double thr_relOffset, thr_absOffset;
	int rr, r;
	int heightPadded, widthPadded, height, width;

	X = mxGetPr(prhs[0]);
	Y = mxGetPr(prhs[1]);
	shield = (bool *) mxGetData(prhs[2]);
	r = mxGetScalar(prhs[3]);
	thr_relOffset = mxGetScalar(prhs[4]);
	thr_absOffset = mxGetScalar(prhs[5]);
	
	heightPadded = mxGetM(prhs[0]);
	widthPadded = mxGetN(prhs[0]);
	
	rr = 2 * r;
	height = heightPadded - rr;
	width = widthPadded - rr;

	plhs[0] = mxCreateDoubleMatrix(height, width, mxREAL);
	double *B = mxGetPr(plhs[0]);
	for (int j = 0; j < width; j++)
	{
		int j_x_height = j*height;
		int jr_x_heightPadded = (j+r)*heightPadded;
		for (int i = 0; i < height; i++)
		{
			int validCount = 0;
			B[j_x_height+i] = 0;
			if (shield[jr_x_heightPadded + i + r])
				continue;

			double xCOffset = X[jr_x_heightPadded + i + r];
			double yCOffset = Y[jr_x_heightPadded + i + r];
			for (int jj = j; jj <= j + rr; jj++)
			{
				int jj_x_heightPadded = jj*heightPadded;
				for (int ii = i; ii <= i + rr; ii++)
				{
					double xOffset = X[jj_x_heightPadded+ii];
					double yOffset = Y[jj_x_heightPadded+ii];
					if (abs(xOffset)+abs(yOffset) >= thr_absOffset)
					{
						validCount++;
						if (abs(xOffset-xCOffset) + abs(yOffset-yCOffset) < thr_relOffset)
						B[j_x_height+i]++;
					}	
				}
			}
			
			if (validCount>(2*r+1)*(2*r+1)*0.2)
				B[j_x_height+i] /= validCount;
			else
				B[j_x_height+i] = 0;
		}
	}
}