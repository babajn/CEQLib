
void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double *b);

void mnewt(int ntrial, double x[], int n, double tolx, double tolf, 
        void (*usrfun)(double *x, int n, double *dvec, double **djac));
