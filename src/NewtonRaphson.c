#include <math.h>
#include "nrutil.h"

#define TINY (1.0e-20)

/************************************************************************
*                          *
*  void ludcmp(a,n,indx,d)                *
*  float **a,*d;                    * 
*  int n,*indx;   n is size of square array        *
*        indx is permutation vector       *
*                          *
*  ludcmp function performs rowwise decomposition        *
*     from numerical recipes in C            *
*                          *
*  ludcmp() takes 4 arguments             *
*  ludcmp() returns no value              *
*                          *
************************************************************************/
// Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition of 
// a rowwise permutation of itself. a and n are input. a is output, arranged as in 
// equation (2.3.14) above; indx[1..n] is an output vector that records the row 
// permutation effected by the partial pivoting; d is output as ±1 depending on 
// whether the number of row interchanges was even or odd, respectively. 
// This routine is used in combination with lubksb to solve linear equations or invert a matrix.
void ludcmp(double **a, int n, int *indx, double *d){
    int imax;
    double big, dum, sum, temp;
    double *vv = dvector(1,n);

    *d = 1.e0;
    for(int i=1; i<=n; i++) {
        big = 0.e0;
        for(int j=1;j<=n;j++)
            if((temp = fabs(a[i][j])) > big) big=temp;
        if(big==0.e0){
           /* nrerror("Singular matrix in routine LUDCMP"); caused EUPDATE event failure in strain preview */
            return ;
        }
        vv[i] = 1.0/big;
    }
    for(int j=1; j<=n; j++) {
        for(int i=1; i<j; i++) {
            sum = a[i][j];
            for(int k=1;k<i;k++)
                sum -= a[i][k]*a[k][j];
            a[i][j] = sum;
        }
        big = 0.e0;
        for(int i=j; i<=n; i++) {
            sum = a[i][j];
            for(int k=1; k<j; k++)
                sum -= a[i][k]*a[k][j];
            a[i][j] = sum;
            if((dum = (vv[i]*fabs(sum))) >= big) {
                big = dum;
                imax = i;
            }
        }
        if(j != imax) {
            for(int k=1; k<=n; k++) {
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            *d = -(*d);
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if(a[j][j] == 0.e0)
            a[j][j] = TINY;
        if(j != n){
            dum = (1.0/(a[j][j]));
            for(int i=j+1; i<=n; i++)
                a[i][j] *= dum;
        }
    }
    free_vector(vv,1,n);
    return(1);
}

/************************************************************************
*									*
*	void lubksb(a,n,indx,b)						*
*	float **a,b[]; 							* 
*	int n,*indx;	n is size of square array			*
*			indx is permutation vector			*
*			returned by ludcmp				*
*									*
*	lubksb function solves a.x = b	from numerical recipes in C	*
*									*
*	lubksb() takes 4 arguments					*
*	lubksb() returns no value					*
*									*
************************************************************************/
// Solves the set of n linear equations A·X = B. Here a[1..n][1..n] is input, 
// not as the matrix A but rather as its LU decomposition, determined by the routine ludcmp. 
// indx[1..n] is input as the permutation vector returned by ludcmp. b[1..n] is input 
// as the right-hand side vector B, and returns with the solution vector X. a, n, and indx 
// are not modified by this routine and can be left in place for successive calls 
// with different right-hand sides b. This routine takes into account the possibility 
// that b will begin with many zero elements, so it is efficient for use in matrix inversion.
void lubksb(double **a, int n, int *indx, double *b){

    int ip, ii = 0;
    double sum;
	for(int i=1; i<=n; i++){
		ip    = indx[i];
		sum   = b[ip];
		b[ip] = b[i];
		if(ii)
			for(int j=ii; j<=i-1; j++)
                sum -= a[i][j]*b[j];
		else if(sum)
            ii = i;
		b[i] = sum;
	}
	for(int i=n; i>=1; i--){
		sum = b[i];
		for(int j=i+1; j<=n; j++)
            sum -= a[i][j]*b[j];
		b[i] = sum/a[i][i];
	}

    return ;
}


// Given an initial guess x[1..n] for a root in n dimensions, 
// take ntrial Newton-Raphson steps to improve the root. 
// Stop if the root converges in either summed absolute variable 
// increments tolx or summed absolute function values tolf.
void mnewt(int ntrial, double x[], int n, double tolx, double tolf, 
        void (*usrfun)(double *x, int n, double *dvec, double **djac)){

    int *indx    = ivector(1,n);
    double *p    = dvector(1,n);
    double *dvec = dvector(1,n);
    double *djac = dmatrix(1,n,1,n);
    double d;

    for(int k=1; k<=ntrial; k++){
        // User function supplies frunction values at x in fvec and Jacobian matrix in fjac.
        usrfun(x,n,dvec,djac);

        // Check function convergence.
        double errf = 0.0;
        for(int i=1;i<=n;i++)
            errf += fabs(dvec[i]);
        if(errf<=tolf){
            //FREERETURN
            free_dmatrix(djac,1,n,1,n);
            free_dvector(dvec,1,n);
            free_dvector(p,1,n);
            free_ivector(indx,1,n);
            return ;
        }

        // Right-hand side of linear equations.
        for(int i=1;i<=n;i++)
            p[i] = -dvec[i];

        // Solve linear equations using LU decomposition.
        ludcmp(djac,n,indx,&d);
        lubksb(djac,n,indx,p);

        // Check root convergence.
        double errx = 0.0;
        // Update solution.
        for(int i=1;i<=n;i++){
            errx += fabs(p[i]);
            x[i] += p[i];
        }
        if(errx<=tolx){
            //FREERETURN
            free_dmatrix(djac,1,n,1,n);
            free_dvector(dvec,1,n);
            free_dvector(p,1,n);
            free_ivector(indx,1,n);
            return ;
        }
    }

    //FREERETURN
    free_dmatrix(djac,1,n,1,n);
    free_dvector(dvec,1,n);
    free_dvector(p,1,n);
    free_ivector(indx,1,n);
    return ;
}
