#include <float.h>
#include <math.h>

#include "ChemicalEquilibriumCooling.h"
#include "CollisionalGasPhaseChemistry.h"
#include "EquilibriumAbundance.h"
#include "FineStructureLines.h"
#include "HIIregionLineCooling.h"
#include "MoleculeChemistry.h"
#include "TemperatureTable.h"

static bool initialized;

static CEQTemperaturePoint InterpolationPoint(int index, double fraction)
{
    /* Four adjacent nodes on a uniform log10(T) grid; store the weights
       once per temperature and evaluate as differences from the lower node. */
    const double throughUpper = fraction*(2.0-fraction);
    const double throughMiddle = fraction*(1.0-fraction);
    return (CEQTemperaturePoint){
        .index = index,
        .before = -throughMiddle*(2.0-fraction)/6.0,
        .upper = throughUpper*(1.0+fraction)/2.0,
        .after = -throughMiddle*(1.0+fraction)/6.0,
    };
}

void CEQFillTemperatureTableWithContext(CEQTemperatureTable *table, double (*evaluate)(double, const void *),
                                        const void *context, const double *breakpoints, size_t count)
{
    const double logSpacing = 1.0/CEQTemperatureBinsPerDecade;
    for(int i=0;i<=CEQTemperatureBins;i++)
        table->sample[i] = evaluate(pow(10.0,i*logSpacing),context);

    /* Use direct evaluations around a fit boundary and wherever a probe
       cannot resolve the coefficient, including an underflow transition. */
    for(int i=1;i<CEQTemperatureBins-1;i++){
        table->interpolate[i] = false;
        const double lower = pow(10.0,(i-1)*logSpacing);
        const double upper = pow(10.0,(i+2)*logSpacing);
        bool usable = true;
        for(size_t j=0;j<count;j++)
            if(breakpoints[j] >= lower && breakpoints[j] <= upper) usable = false;
        for(int j=i-1;j<=i+2;j++){
            const double value = table->sample[j];
            if(!isfinite(value) || (value == 0.0) != (table->sample[i] == 0.0) ||
               (value != 0.0 && fabs(value) < DBL_MIN)) usable = false;
        }
        if(!usable) continue;
        for(int j=1;j<=3;j++){
            const double fraction = 0.25*j;
            const double exact = evaluate(pow(10.0,(i+fraction)*logSpacing),context);
            const CEQTemperaturePoint point = InterpolationPoint(i,fraction);
            const double estimated = CEQInterpolateTemperature(table->sample+i,&point);
            if(!isfinite(exact) || !isfinite(estimated) ||
               (fabs(exact) < DBL_MIN ? estimated != exact : fabs(estimated-exact) > 5.0e-9*fabs(exact))){
                usable = false;
                break;
            }
        }
        table->interpolate[i] = usable;
    }
}

static double EvaluateWithoutContext(double temperature, const void *context)
{
    double (*const *evaluate)(double) = context;
    return (*evaluate)(temperature);
}

void CEQFillTemperatureTable(CEQTemperatureTable *table, double (*evaluate)(double),
                             const double *breakpoints, size_t count)
{
    CEQFillTemperatureTableWithContext(table,EvaluateWithoutContext,&evaluate,breakpoints,count);
}

void CEQInitializeNonEquilibriumTemperatureTables(void)
{
    if(initialized) return;
    CEQInitializeChemistryTemperatureTables();
    CEQInitializeAtomicCoolingTemperatureTables();
    CEQInitializeMolecularCoolingTemperatureTables();
    CEQInitializeOxygenCoolingTemperatureTables();
    CEQInitializeFineStructureTemperatureTables();
    initialized = true;
}

bool CEQLocateLogTemperature(double logTemperature, CEQTemperaturePoint *point)
{
    if(!initialized || !EQCLibRunParameters.UseNonEquilibriumTemperatureTables) return false;
    const double scaled = logTemperature*CEQTemperatureBinsPerDecade;
    if(!(scaled >= 1.0 && scaled < CEQTemperatureBins-1)) return false;
    const int index = (int)scaled;
    *point = InterpolationPoint(index,scaled-index);
    return true;
}

bool CEQLocateTemperature(double temperature, CEQTemperaturePoint *point)
{
    return temperature > 1.0 && CEQLocateLogTemperature(log10(temperature),point);
}
