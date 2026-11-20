#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "Astro.h"
#include "Constants.h"
#include "ChemicalEquilibriumCooling.h"
#include "ChemicalElements.h"
#include "CollisionalGasPhaseChemistry.h"
#include "PhotochemicalGasPhaseChemistry.h"
#include "GrainSurfaceChemistry.h"
#include "MoleculeChemistry.h"
#include "CosmicRay.h"
#include "Xray.h"
    
#define TURBULENCE_VELOCITY (3.0)   // km/s
double b_turb = TURBULENCE_VELOCITY;    // km/s


void EquilibriumHydrogenFractionalAbundance(StructCEQData *PhydroPP){

    double ne     = PhydroPP->nel;
    double nH     = PhydroPP->n[ipHYDROGEN];
	double nHII   = PhydroPP->nion [ipHYDROGEN];
    double nHI    = PhydroPP->natom[ipHYDROGEN];
    double nH2    = PhydroPP->nmol [ipH2];
    double T      = PhydroPP->T;
    double NLyn   = PhydroPP->NLyn;
    double G0     = PhydroPP->G0;
    double X0     = PhydroPP->X0;
    double Ntot   = PhydroPP->Nt;
    double Nh2    = Ntot*(nH2/nH);
    double Zd     = PhydroPP->Zd;
    double xe     = ne/nH;
    double Tdust  = ApproximateDustTemperature(Ntot,G0);
    //double b_turb = TURBULENCE_VELOCITY;    // km/s
    double b5     = sqrt(SQ(12.90*sqrt(1.e-4*T))+SQ(b_turb));    // eq(6.34) of Draine's textbook

    double dnHI_phot  = HIPhotoIonizationRate(NLyn);
    double dnHI_coli  = ne*HICollisionalIonizationRate(T);
    double dnHI_cr    = X0*CosmicRayIonizationRateHI(xe);
    double dnHI_xr    = X0*SoftXrayIonizationRateHI(xe,Ntot);
    double dnHII_rec  = ne*HIIRecombinationRate(T);
    double dnHII_dust = nH*HIIRecombinationRateOnDust(T,G0,ne,Zd); 
    double HIItoHI    = (dnHI_phot + dnHI_coli + dnHI_cr + dnHI_xr)/(dnHII_rec + dnHII_dust);
    double fHII, fHI, fH2;
    if(CEQLibRunParameters.UseH2Formation==true){
        double dnH2_form  = nH*H2FormationRate(T,Zd,Tdust);
        double dnH2_diss  = H2PhotoDissociationRate(G0,Ntot,Nh2,b5,Zd);
        double dnH2_dest  = H2CollisionalDissociationRate(T,nH,nHI,nH2,ne);
        double H2toHI     = dnH2_form/(dnH2_diss + dnH2_dest);
        double HIItoHt    = MIN(1.0,HIItoHI/(HIItoHI+2*H2toHI+1.0));
        double HItoHt     = MIN(1.0,1.0    /(HIItoHI+2*H2toHI+1.0));
        fHII = HIItoHt;
        fHI  = HItoHt;
        fH2  = 0.5*MAX(0.0,1.e0 - fHII - fHI);
    }else{
        double HIItoHt    = MIN(1.0,HIItoHI/(HIItoHI+1.0));
        double HItoHt     = MIN(1.0,1.0    /(HIItoHI+1.0));
        fHII = HIItoHt;
        fHI  = MAX(0.0,1.e0 - fHII);
        fH2  = 0.0;
    }
    if(fHII<0||fHI<0||fH2<0){
        fprintf(stderr,"#%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
        fprintf(stderr,"  fHII=%e fHI=%e fH2=%e\n",fHII,fHI,fH2);
        exit(1);
    }
	PhydroPP->nion[ipHYDROGEN]  = nH*fHII;
    PhydroPP->natom[ipHYDROGEN] = nH*fHI ;
    PhydroPP->nmol[ipH2]        = nH*fH2 ;
    return ;
}

//void EquilibriumHeliumFractionalAbundance(double T, double nH, double ne, double G0, double Ntot, double Z, double *fHeII, double *fHeI){
//    double xe          = ne/nH;
//    double X0          = G0;
//    double dnHeI_coli  = ne*HeICollisionalIonizationRate(T);
//    double dnHeI_cr    = X0*CosmicRayIonizationRateHeI(xe);
//    double dnHeI_xr    = X0*SoftXrayIonizationRateHeI(xe,Ntot);
//    double dnHeII_rec  = ne*(HeIIRecombinationRateCaseA(T) + HeIIDielectronicRecombinationRate(T));
//    double dnHeII_dust = nH*HeIIRecombinationRateOnDust(T,G0,ne,Z); 
//
//    double HeIItoHeI   = (dnHeI_coli + dnHeI_cr + dnHeI_xr)/(dnHeII_rec + dnHeII_dust);
//    *fHeII = HeIItoHeI/(1.0+HeIItoHeI);
//    *fHeI  = MAX(0.0,1.0 - (*fHeII));
//    if(*fHeII<0||*fHeI<0){
//        fprintf(stderr,"#%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
//        fprintf(stderr,"  fHeII=%e fHeI=%e\n",*fHeII,*fHeI);
//        exit(1);
//    }
//    
//    return ;
//}

void EquilibriumHeliumFractionalAbundance(StructCEQData *PhydroPP){
    double ne     = PhydroPP->nel;
    double nH     = PhydroPP->n[ipHYDROGEN];
    double nHe    = PhydroPP->n[ipHELIUM];
	double nHeIII = PhydroPP->nHeIII;
	double nHeII  = PhydroPP->nion [ipHELIUM];
    double nHeI   = PhydroPP->natom[ipHELIUM];
    double T      = PhydroPP->T;
    double NLyn   = PhydroPP->NLyn; 
    double G0     = PhydroPP->G0;
    double X0     = PhydroPP->X0;
    double Ntot   = PhydroPP->Nt;
    double Zd     = PhydroPP->Zd;
    double xe     = ne/nH;
    
    double dnHeI_phti   = HeIPhotoIonizationRate(NLyn);
    double dnHeI_coli   = ne*HeICollisionalIonizationRate(T);
    double dnHeI_cr     = X0*CosmicRayIonizationRateHeI(xe);
    double dnHeI_xr     = X0*SoftXrayIonizationRateHeI(xe,Ntot);
    double dnHeII_rec   = ne*(HeIIRecombinationRate(T) + HeIIDielectronicRecombinationRate(T));
    double dnHeII_dust  = nH*HeIIRecombinationRateOnDust(T,G0,ne,Zd); 

    double dnHeII_phti  = HeIIPhotoIonizationRate(NLyn);
    double dnHeII_coli  = ne*HeIICollisionalIonizationRate(T);
    double dnHeII_cr    = 0.e0; // OK
    double dnHeII_xr    = 0.e0; // 0K
    double dnHeIII_rec  = ne*HeIIIRecombinationRate(T);
    double dnHeIII_dust = 0.e0; // OK
    
    double HeIItoHeI   = (dnHeI_phti  + dnHeI_coli  + dnHeI_cr  + dnHeI_xr )/(dnHeII_rec  + dnHeII_dust );
    double HeIIItoHeII = (dnHeII_phti + dnHeII_coli + dnHeII_cr + dnHeII_xr)/(dnHeIII_rec + dnHeIII_dust);
    double fHeIII = HeIIItoHeII/(HeIIItoHeII + 1.e0 + 1.e0/HeIItoHeI);
    double fHeII  =        1.e0/(HeIIItoHeII + 1.e0 + 1.e0/HeIItoHeI);
    double fHeI   = MAX(0.0,1.e0 - fHeII - fHeIII);
    if(fHeIII<0||fHeII<0||fHeI<0){
        fprintf(stderr,"#%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
        //fprintf(stderr,"  ne=%e nH=%e T=%e\n",ne,nH,T);
        //fprintf(stderr,"  G0=%e X0=%e NLyn=%e\n",G0,X0,PhydroPP->NLyn);
        //fprintf(stderr,"  Nt=%e Z=%e\n",Ntot,Z);
        //fprintf(stderr,"  dnHeI_phti=%e dnHeI_coli=%e dnHeI_cr=%e dnHeI_xr=%e\n",dnHeI_phti,dnHeI_coli,dnHeI_cr,dnHeI_xr);
        fprintf(stderr,"  fHeIII=%e fHeII=%e fHeI=%e\n",fHeIII,fHeII,fHeI);
        exit(1);
    }
    PhydroPP->nHeIII          = nHe*fHeIII;
    PhydroPP->nion[ipHELIUM]  = nHe*fHeII ;
    PhydroPP->natom[ipHELIUM] = nHe*fHeI  ;

    return ;
}


//// ionization state of oxygen (charge transfer)
//// ref. eq(2.38) in Osterbrock & Ferland (2006)
////      eq(14.32) in Draine (2011)
#define dE_OI_HI (0.19*EV_TO_ERG)
void EquilibriumOxygenFractionalAbundanceChargeTransfer(StructCEQData *PhydroPP){
    double kT   = BOLTZMANN_CONSTANT_CGS*PhydroPP->T;
    double nHI  = PhydroPP->nion [ipHYDROGEN];
    double nHII = PhydroPP->natom[ipHYDROGEN];
    double nO   = PhydroPP->n[ipOXYGEN];
    double nCO  = PhydroPP->nmol [ipCO];
    double RateRatio = 9.0*exp(dE_OI_HI/kT)/8.0;
    double fOIII = 0.e0;
    double fOII  = nHII/(nHII + RateRatio*nHI);
    double fOI   = MAX(0.0,1.0 - fOII);
    if(fOII<0||fOI<0){
        fprintf(stderr,"#%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
        fprintf(stderr,"fOII = %e fOI = %e T = %e nHI = %e nHII = %e\n",
                fOII,fOI,kT/BOLTZMANN_CONSTANT_CGS,nHI,nHII);
        exit(1);
    }
    PhydroPP->nOIII           = MAX(0.0,nO-nCO)*fOIII;
    PhydroPP->nion[ipOXYGEN]  = MAX(0.0,nO-nCO)*fOII ;
    PhydroPP->natom[ipOXYGEN] = MAX(0.0,nO-nCO)*fOI  ;
    return ;
}

void EquilibriumOxygenFractionalAbundance(StructCEQData *PhydroPP){

    double ne     = PhydroPP->nel;
    double nH     = PhydroPP->n[ipHYDROGEN];
    double nHI    = PhydroPP->nion [ipHYDROGEN];
    double nHII   = PhydroPP->natom[ipHYDROGEN];
    double nH2    = PhydroPP->nmol [ipH2];
    double nO     = PhydroPP->n[ipOXYGEN];
    double nOII   = PhydroPP->nion [ipOXYGEN];
    double nOI    = PhydroPP->natom[ipOXYGEN];
    double nCO    = PhydroPP->nmol [ipCO];
    double T      = PhydroPP->T;
    double NLyn   = PhydroPP->NLyn;
    double G0     = PhydroPP->G0;
    double X0     = PhydroPP->X0;
    double Ntot   = PhydroPP->Nt;
    double Zd     = PhydroPP->Zd;
    double xe     = ne/nH;
    
    double dnOI_phti   = OIPhotoionizationRate(NLyn);
    double dnOI_coli   = ne*OICollisionalIonizationRate(T);
    double dnOI_cr     = X0*CosmicRayIonizationRateOI(xe);
    double dnOI_xr     = 0.e0;
    double dnOII_rec   = ne*OIIRecombinationRate(T);
    double dnOII_dust  = nH*OIIRecombinationRateOnDust(T,G0,ne,Zd);
    double dnOI_chtr   = 0.e0; //nHII*OIChargeTransferRateToHII(T);
    double dnOII_chtr  = 0.e0; //nHI *OIIChargeTransferRateToHI(T);

    double dnOII_phti  = OIIPhotoionizationRate(NLyn);
    double dnOII_coli  = ne*OIICollisionalIonizationRate(T);
    double dnOII_cr    = 0.e0;
    double dnOII_xr    = 0.e0;
    double dnOIII_rec  = ne*OIIIRecombinationRate(T);
    double dnOIII_dust = 0.e0;
    
    double OIItoOI   = (dnOI_phti  + dnOI_coli  + dnOI_cr  + dnOI_xr  + dnOI_chtr)/(dnOII_rec  + dnOII_dust + dnOII_chtr);
    double OIIItoOII = (dnOII_phti + dnOII_coli + dnOII_cr + dnOII_xr)/(dnOIII_rec + dnOIII_dust);
    double fOIII = OIIItoOII/(OIIItoOII + 1.e0 + 1.e0/OIItoOI);
    double fOII  =      1.e0/(OIIItoOII + 1.e0 + 1.e0/OIItoOI);
    double fOI   = MAX(0.0,1.e0 - fOII - fOIII);
    if(fOIII<0||fOII<0||fOI<0){
        fprintf(stderr,"#%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
        fprintf(stderr,"fOIII = %e fOII = %e fOI = %e \n",fOIII,fOII,fOI);
        exit(1);
    }
    PhydroPP->nOIII           = MAX(0.0,nO-nCO)*fOIII;
    PhydroPP->nion[ipOXYGEN]  = MAX(0.0,nO-nCO)*fOII ;
    PhydroPP->natom[ipOXYGEN] = MAX(0.0,nO-nCO)*fOI  ;

    return ;
}


//#ifdef CHEMISTRY_CO_KETO_CASELLI_2008
// equilibrium carbon chemistry
// ref. Keto & Caselli, ApJ, 683, 238, (2008)
//      see footnote 2 of Glover & Clark, MNRAS, 421, 116, (2012)
static double COFormationRateCH2(double nH2, double T){
    // CH2 formation rate coefficient from CII and H2 (via radiative association to form CH2+)
    //double kRA = 6.e-16; // ref.  Table 12 of Tielens & Hollenbach (1985)
    double Cp = CEQLibRunParameters.CpForCO;
    double kRA = 4.e-16*pow(T/300.0,-0.2);  // k148 in Table B1 of Glover et al. (2010)
    return kRA*Cp*nH2; // [/s] 
}

static double COtoCIIratio(double nH2, double T, double G0, double Ntot, double N2, double Nco, double Z){
    // eq.(14) in Glover & Clark, MNRAS, 421, 116, (2012)
    //double b_turb = TURBULENCE_VELOCITY;    // km/s
    double b5   = sqrt(SQ(12.90*sqrt(1.e-4*T))+SQ(b_turb));    // eq(6.34) of Draine's textbook
    double r_CO = COPhotoDissociationRate(G0,Ntot,Z)*COShieldingFactor(N2,Nco,b5);
    return COFormationRateCH2(nH2,T)/r_CO;
}

static double CItoCIIratio(double nH2, double T, double G0, double Ntot, double Z){
    // eq.(15) in Glover & Clark, MNRAS, 421, 116, (2012)
    double r_C = 3.1e-10*G0*DustShieldingFactor(Ntot,Z,3.0); // R171 in Table B2 of Glover et al. 2010
    return COFormationRateCH2(nH2,T)/r_C;
}

// eq(11)-(13) in Glover & Clark, MNRAS, 421, 116, (2012)
void EquilibriumCarbonFractionalAbundance(StructCEQData *PhydroPP){
    double ne   = PhydroPP->nel;
    double nH   = PhydroPP->n[ipHYDROGEN];
    double nH2  = PhydroPP->nmol [ipH2];
    double nC   = PhydroPP->n[ipCARBON];
    double nCII = PhydroPP->nion [ipCARBON];
    double nCI  = PhydroPP->natom[ipCARBON];
    double T    = PhydroPP->T;
    double G0   = PhydroPP->G0;
    double X0   = PhydroPP->X0;
    double Ntot = PhydroPP->Nt;
    double Zd   = PhydroPP->Zd;
    double xe   = ne/nH;
    double fCIII = 0.e0, fCII = 0.e0, fCI = 0.e0, fCO = 0.e0;
    if(CEQLibRunParameters.UseCOFormation==true){
        double nCO = PhydroPP->nmol [ipCO];
        double Nh2 = Ntot*(nH2/nH);
        double Nco = Ntot*nCO/nH;
        double COtoCII = COtoCIIratio(nH2,T,G0,Ntot,Nh2,Nco,Zd); // eq.(14)
        double CItoCII = CItoCIIratio(nH2,T,G0,Ntot,Zd);         // eq.(15)
        double CIItoCtotal = 1.e0/(1.0 + COtoCII + CItoCII);    // eq.(11)
        fCIII = 0.e0;
        fCII  = CIItoCtotal;
        fCI   = CItoCII*CIItoCtotal;
        fCO   = MAX(0.0,1.e0 - fCII - fCI);
        if(fCIII<0||fCII<0||fCI<0||fCO<0){
            //double b_turb = 0.e0;    // km/s
            double b5     = sqrt(SQ(12.90*sqrt(1.e-4*T))+SQ(b_turb));    // eq(6.34) of Draine's textbook
            fprintf(stderr,"#%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
            fprintf(stderr,"  COtoCII=%e CItoCII=%e \n",COtoCII,CItoCII);
            fprintf(stderr,"  nH2=%e T=%e G0=%e Ntot=%e Nh2=%e Nco=%e Zd=%e \n",nH2,T,G0,Ntot,Nh2,Nco,Zd);
            fprintf(stderr,"  COSelfShielding=%e b5=%e\n",COShieldingFactor(Nh2,Nco,b5),b5);
            fprintf(stderr,"  fCIII=%e fCII=%e fCI=%e fCO=%e\n",fCIII,fCII,fCI,fCO);
            exit(1);
        }
    //} else{
    //} if(CEQLibRunParameters.UseCOFormation==false||fCO<1.e-9){
    } if(CEQLibRunParameters.UseCOFormation==false||T>1.e+4){
        double NLyn       = PhydroPP->NLyn;
        double xe         = ne/nH;
        //double Av         = 0.0; //VisualExtinction(Ntot,Z);    for Wolfire et al. 1995
        double Av         = VisualExtinction(Ntot,Zd);
        
        double dnCI_phti  = 0.e0;
        if(NLyn==0.e0)  dnCI_phti = CIFUVPhotoionizationRate(G0,Av);
        else            dnCI_phti = CIPhotoionizationRate(NLyn);
        double dnCI_coli  = ne*CICollisionalIonizationRate(T);
        double dnCI_cr    = X0*CosmicRayIonizationRateCI(xe);
        double dnCI_xr    = X0*SoftXrayIonizationRateCI(xe,Ntot);
        double dnCII_rec  = ne*CIIRecombinationRate(T);
        double dnCII_dust = nH*CIIRecombinationRateOnDust(T,G0,ne,Zd); 

        double dnCII_phti  = CIIPhotoionizationRate(NLyn);
        double dnCII_coli  = ne*CIICollisionalIonizationRate(T);
        double dnCII_cr    = 0.e0;
        double dnCII_xr    = 0.e0;
        double dnCIII_rec  = ne*CIIIRecombinationRate(T);
        double dnCIII_dust = 0.e0;
       
        double CIItoCI   = (dnCI_phti  + dnCI_coli  + dnCI_cr  + dnCI_xr )/(dnCII_rec  + dnCII_dust );
        double CIIItoCII = (dnCII_phti + dnCII_coli + dnCII_cr + dnCII_xr)/(dnCIII_rec + dnCIII_dust);

        //fCII = CIItoCI/(1.0+CIItoCI);
        //fCI  = MAX(0.0,1.0 - fCII);
        //fCO  = 0.e0;
        fCIII = CIIItoCII/(CIIItoCII + 1.e0 + 1.e0/CIItoCI);
        fCII  =      1.e0/(CIIItoCII + 1.e0 + 1.e0/CIItoCI);
        fCI   = MAX(0.0,1.e0 - fCII - fCIII);
        fCO   = 0.e0;
        if(fCIII<0||fCII<0||fCI<0||fCO<0){
            fprintf(stderr,"#%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
            fprintf(stderr,"  fCIII=%e fCII=%e fCI=%e fCO=%e\n",fCIII,fCII,fCI,fCO);
            exit(1);
        }
    }
	PhydroPP->nCIII           = nC*fCIII;
	PhydroPP->nion[ipCARBON]  = nC*fCII ;
    PhydroPP->natom[ipCARBON] = nC*fCI  ;
    PhydroPP->nmol[ipCO]      = nC*fCO  ;
    return ;
}

void EquilibriumSiliconFractionalAbundance(StructCEQData *PhydroPP){
    double ne     = PhydroPP->nel;
    double nH     = PhydroPP->n[ipHYDROGEN];
    double nSi    = PhydroPP->n[ipSILICON];
    double nSiII  = PhydroPP->nion [ipSILICON];
    double nSiI   = PhydroPP->natom[ipSILICON];
    double T      = PhydroPP->T;
    double NLyn   = PhydroPP->NLyn;
    double G0     = PhydroPP->G0;
    double X0     = PhydroPP->X0;
    double Ntot   = PhydroPP->Nt;
    double Zd     = PhydroPP->Zd;
    double xe     = ne/nH;

    double dnSiI_phti   = SiIPhotoionizationRate(NLyn);
    double dnSiI_coli   = ne*SiICollisionalIonizationRate(T);
    double dnSiI_cr     = X0*CosmicRayIonizationRateSiI(xe);
    double dnSiI_xr     = X0*SoftXrayIonizationRateSiI(xe,Ntot);
    double dnSiII_rec   = ne*SiIIRecombinationRate(T);
    double dnSiII_dust  = nH*SiIIRecombinationRateOnDust(T,G0,ne,Zd); 

    double dnSiII_phti  = SiIIPhotoionizationRate(NLyn);
    double dnSiII_coli  = ne*SiIICollisionalIonizationRate(T);
    double dnSiII_cr    = X0*CosmicRayIonizationRateSiII(xe);
    double dnSiII_xr    = X0*SoftXrayIonizationRateSiII(xe,Ntot);
    double dnSiIII_rec  = ne*SiIIIRecombinationRate(T);
    double dnSiIII_dust = 0.e0; // OK
    
    double SiIItoSiI   = (dnSiI_phti  + dnSiI_coli  + dnSiI_cr  + dnSiI_xr )/(dnSiII_rec  + dnSiII_dust );
    double SiIIItoSiII = (dnSiII_phti + dnSiII_coli + dnSiII_cr + dnSiII_xr)/(dnSiIII_rec + dnSiIII_dust);
   
    double fSiIII = SiIIItoSiII/(SiIIItoSiII + 1.e0 + 1.e0/SiIItoSiI);
    double fSiII  =        1.e0/(SiIIItoSiII + 1.e0 + 1.e0/SiIItoSiI);
    double fSiI   = MAX(0.0,1.e0 - fSiII - fSiIII);
    if(fSiIII<0||fSiII<0||fSiI<0){
        fprintf(stderr,"#%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
        fprintf(stderr,"  fSiIII=%e fSiII=%e fSiI=%e\n",fSiIII,fSiII,fSiI);
        exit(1);
    }
    PhydroPP->nSiIII            = nSi*fSiIII;
    PhydroPP->nion[ipSILICON]   = nSi*fSiII ;
    PhydroPP->natom[ipSILICON]  = nSi*fSiI  ;
    return ;
}

