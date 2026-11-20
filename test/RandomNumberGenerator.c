#include <math.h>
#include <gsl/gsl_rng.h>
#include "RandomNumberGenerator.h"

/* Random Generator */
gsl_rng *RandomGenerator;
/* Random Generator */

void InitializeRandomGenerator(unsigned long int RandomSeed){

    AllocateRandomGenerator();
    //RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(RandomGenerator,RandomSeed);

    fprintf(stderr,"Random generator name : %s\n",gsl_rng_name(RandomGenerator));
    fprintf(stderr," Seed %ld\n",RandomSeed);
    fprintf(stderr," Max and min value of the random numbers: %ld <-> %ld\n",
            gsl_rng_min(RandomGenerator),gsl_rng_max(RandomGenerator));


    return;
}

void AllocateRandomGenerator(void){
    RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
    return;
}

/***** Gaussian with mean=0.0 and dispersion=1.0 by Box-Muller method *****/
double gaussian(void){
    double x, y, z, r2;

    do{
        x = 2.0*gsl_rng_uniform(RandomGenerator) - 1.0;
        y = 2.0*gsl_rng_uniform(RandomGenerator) - 1.0;
        r2 = x*x + y*y;
    }while(r2 >= 1.0 || r2 == 0.0);
    z = sqrt(-2.0*log(r2)/r2)*x;
    
    return (z);
}


