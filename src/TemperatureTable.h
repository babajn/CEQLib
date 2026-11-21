#ifndef CEQ_TEMPERATURE_TABLE_H
#define CEQ_TEMPERATURE_TABLE_H

#include <math.h>
#include <stdbool.h>
#include <stddef.h>

enum {
    CEQTemperatureBinsPerDecade = 1024,
    CEQTemperatureBins = 8*CEQTemperatureBinsPerDecade,
};

typedef struct {
    double sample[CEQTemperatureBins+1];
    bool interpolate[CEQTemperatureBins];
} CEQTemperatureTable;

typedef struct {
    int index;
    double before, upper, after;
} CEQTemperaturePoint;

void CEQFillTemperatureTable(CEQTemperatureTable *table, double (*evaluate)(double),
                             const double *breakpoints, size_t count);
void CEQFillTemperatureTableWithContext(CEQTemperatureTable *table, double (*evaluate)(double, const void *),
                                        const void *context, const double *breakpoints, size_t count);
void CEQInitializeNonEquilibriumTemperatureTables(void);
bool CEQLocateTemperature(double temperature, CEQTemperaturePoint *point);
bool CEQLocateLogTemperature(double logTemperature, CEQTemperaturePoint *point);

static inline double CEQInterpolateTemperature(const double *samples, const CEQTemperaturePoint *point)
{
    const double center = samples[0];
    return center + (point->before*(samples[-1]-center) + point->upper*(samples[1]-center)
                     + point->after*(samples[2]-center));
}

/* An unusable interval evaluates the original rate. Pass NULL for direct to
   handle fallback at the call site, for instance when it supplies a vector. */
static inline double CEQReadTemperatureTable(const CEQTemperatureTable *table, const CEQTemperaturePoint *point,
                                              double temperature, double (*direct)(double))
{
    if(point != NULL && table->interpolate[point->index])
        return CEQInterpolateTemperature(table->sample+point->index,point);
    return direct != NULL ? direct(temperature) : NAN;
}

#endif
