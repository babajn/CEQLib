#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Astro.h"
#include "ChemicalElements.h"
#include "Bremsstrahlung.h"

//// Bremsstrahlung Cooling [erg/cc/s]
// ref. eq.(16) of Cen (1992) 
#define GauntFactor (1.5)
double CoolingBremsstrahalung(StructCEQData *PhydroPP){
    return 1.42e-27*GauntFactor*sqrt(PhydroPP->T)*PhydroPP->nel
            *(PhydroPP->nion[ipHYDROGEN]+PhydroPP->nion[ipHELIUM]+4.0*PhydroPP->nHeIII); 
}
