#ifndef __EQCDATASTRUCTURE_H_INCLUDED__
#define __EQCDATASTRUCTURE_H_INCLUDED__
typedef struct StructEQCData_tag{
    int GlobalID;
    double t;
    double dt;
    double Mass;
    double Pos[3];
    double Vel[3];
    double Rho;
    double n[6];
    double nHeIII;
    double nCIII;
    double nOIII;
    double nSiIII;
    double nion[6];
    double natom[6];
    double nmol[2];
    double nel;
    double U;          // Inner energy of gas.
    double Du;
    double DQheat;
    double T;
    double invT;
    double logT;
    double Td;          // Positive dust temperature in K, otherwise use the library estimate.
    double Z;           // Metallicity of gas.
    double Zd;          // Dust/Gas mass ratio
    double NLyn;        // Lymann photons (>13.6 eV) flux [/s/cm^2]
    double G0;          // FUV
    double X0;          // XR & CR
    double Nt;          // [H/cm^2]
    //double Kernel;      // Size of kernel. 0<- rij/Kernel ->2
    double MeanMolecularWeight;
    double DegreeOfFreedom;
    double Gamma;
    double Gm1;
    double GGm1;

    // heating rate [erg/s/cc]
    double H_tot;
    double H_PI;
    double H_PE;
    double H_CR;    
    double H_XR;    
    double H_H2;    

    // cooling rate, i.e. angular-integrated emission coefficient (= 4*pi*j) [erg/s/cc]
    double L_tot;
    double L_res;    
    double L_ff;    
    double L_HI;    
    double L_HeI;    
    double L_HeII;    
    double L_CII158;  
    double L_CI610;
    double L_CI370;
    double L_OIII;      // [O III] line cooling included in L_tot, with H II weighting
    double L_OII;       // [O II] line cooling included in L_tot, with H II weighting
    double L_OI63;      // [O I]63um 
    double L_OI146;     // [O I]146um
    double L_SiII35;
    double L_DustRec;
    double L_H2;
    double L_CO;
    
    // Halpha emission
    double L_Halpha;    

    // Hi 21cm
    // 4pi * line-integrated emissivity [erg/s/cc]
    double L_HI21cm;    
    // absorption coefficient [1/cm]
    double k_HI21cm;    
} StructEQCData;
#endif
