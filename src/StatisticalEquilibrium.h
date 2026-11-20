double EinsteinCoefficientB10(double A10, double E10);
double EinsteinCoefficientB01(double B10, double g1g0);

void SolveTwoLevelSystemStatisticalEquilibrium(
        double E10, double g1g0,
        double I10,double A10, double B10, double B01, 
        int NcolSpecies, double ncol[], double q10[], 
        double n, double kT, double *n0, double *n1);
void SolveThreeLevelSystemStatisticalEquilibrium(
        double E10, double E20, double E21,
        double g1g0,double g2g0,double g2g1,
        double I10, double I20, double I21,
        double A10, double A20, double A21,
        double B10, double B20, double B21,
        double B01, double B02, double B12,
        int NcolSpecies, double ncol[], double q10[], double q20[], double q21[],
        double n, double kT, double *n0, double *n1, double *n2);
