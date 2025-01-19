/*
  Optimized double class functions for MATLAB on x86 computers.
  Copyright � Marcel Leutenegger, 2003-2008, �cole Polytechnique F�d�rale de Lausanne (EPFL),
  Laboratoire d'Optique Biom�dicale (LOB), BM - Station 17, 1015 Lausanne, Switzerland.

      This library is free software; you can redistribute it and/or modify it under
      the terms of the GNU Lesser General Public License as published by the Free
      Software Foundation; version 2.1 of the License.

      This library is distributed in the hope that it will be useful, but WITHOUT ANY
      WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
      PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

      You should have received a copy of the GNU Lesser General Public License along
      with this library; if not, write to the Free Software Foundation, Inc.,
      51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

/****************************************************************************\
*									     *
*	      Construct complex data from real and imaginary parts	     *
*									     *
*									     *
*			Marcel Leutenegger � 12.5.2008			     *
*									     *
\****************************************************************************/

#include "mex.h"
#define	O	plhs[0]
#define	S	prhs[0]
#define	T	prhs[1]

bool fdim(int m, int n, const int* d, const int* e);
void fmov(double* or, const double* sr, int n);
void fset(double* or, const double* sr, int n);


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{	double* or;
	int m, n;
	if (nlhs > 1) mexErrMsgTxt("Too many output arguments.");
	switch(nrhs)
	{default:
		mexErrMsgTxt("Wrong number of input arguments.");
	 case 0:
		mexPrintf("\nConstruct complex data from real and imaginary parts.\n\n\tMarcel Leutenegger � 12.5.2008\n\n");
		break;
	 case 1:
		if (mxIsComplex(S)) mexErrMsgTxt("Arguments must be real.");
		O=mxCreateDoubleMatrix(0,0,mxREAL);
		mxSetDimensions(O,mxGetDimensions(S),mxGetNumberOfDimensions(S));
		n=mxGetNumberOfElements(S);
		if (n > 0)
		{	mxSetPi(O,mxCalloc(n,sizeof(double)));
			mxSetPr(O,or=mxMalloc(n*sizeof(double)));
			fmov(or,mxGetPr(S),n);
		}
		break;
	 case 2:
		if (mxIsComplex(S) || mxIsComplex(T)) mexErrMsgTxt("Arguments must be real.");
		O=mxCreateDoubleMatrix(0,0,mxREAL);
		m=mxGetNumberOfElements(S);
		n=mxGetNumberOfElements(T);
		if (m > 0 && n > 0)
		{	if (m == 1)	// scalar,matrix
			{	mxSetDimensions(O,mxGetDimensions(T),mxGetNumberOfDimensions(T));
				mxSetPr(O,or=mxMalloc(n*sizeof(double)));
				fset(or,mxGetPr(S),n);
				mxSetPi(O,or=mxMalloc(n*sizeof(double)));
				fmov(or,mxGetPr(T),n);
				break;
			}
			if (n == 1)	// matrix,scalar
			{	mxSetDimensions(O,mxGetDimensions(S),mxGetNumberOfDimensions(S));
				mxSetPr(O,or=mxMalloc(m*sizeof(double)));
				fmov(or,mxGetPr(S),m);
				mxSetPi(O,or=mxMalloc(m*sizeof(double)));
				fset(or,mxGetPr(T),m);
				break;
			}
			if (m == n)	// matrix,matrix
			{	const int* d=mxGetDimensions(S);
				m=mxGetNumberOfDimensions(S);
				if (fdim(m,mxGetNumberOfDimensions(T),d,mxGetDimensions(T)))
				{	mxSetDimensions(O,d,m);
					mxSetPr(O,or=mxMalloc(n*sizeof(double)));
					fmov(or,mxGetPr(S),n);
					mxSetPi(O,or=mxMalloc(n*sizeof(double)));
					fmov(or,mxGetPr(T),n);
					break;
				}
			}
			mexErrMsgTxt("Incompatible dimensions.");
		}
	}
}
