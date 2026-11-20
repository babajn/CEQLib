/* program for reading */

#include <stdio.h>
#include <stdlib.h>

main()
{
	char inputfile[80], scenarios[200];	
	int  j, k, h;
	int  nbtimesteps=10,nbcontinuum=10, nblines=10;
	double wave_conti[1500], wave_lines[100], flux_conti[1500], flux_lines[100];
	double age, a, M_galaxy, M_star, M_WD, M_NS, M_substar, M_gas, metal_ISM, metal_star_mass, metal_star_BL;
	double BL, tau_V, dust, SFR, Ly, SN2, SN1a, age_star_mass, age_star_BL;
	FILE *FRV, *synthspecflux[10], *scenariofile;
	
	FRV=fopen("specE.dat","r");
	scenariofile=fopen("out_scenario.dat","w");
	synthspecflux[1]=fopen("out_spec.dat","w");
	/*Evolutionary Scenarioの読み込み*/
	printf("ok\n");
	for(j=1;j<=31;j++){
		fgets(scenarios, 200, FRV);
		fputs(scenarios, scenariofile);
	}
	printf("ok\n");
	/*timestep、波長サンプル数の読み込み*/
	fscanf(FRV,"%d %d %d",&nbtimesteps,&nbcontinuum,&nblines);
	/*波長の読み込み*/
	for(j=1; j<=nbcontinuum; j++)
			fscanf(FRV,"%lf",&wave_conti[j]);
	for(j=1; j<=nblines; j++)
			fscanf(FRV,"%lf",&wave_lines[j]);
	/*nbtimesteps回分読み込む*/
	a=0;
	for(j=1; j<=nbtimesteps; j++)
	{
		/*propertyの読み込み*/
		fscanf(FRV,"%lf",&age);
		fscanf(FRV,"%lf",&M_galaxy);
		fscanf(FRV,"%lf",&M_star);
		fscanf(FRV,"%lf",&M_WD);
		fscanf(FRV,"%lf",&M_NS);
		fscanf(FRV,"%lf",&M_substar);
		fscanf(FRV,"%lf",&M_gas);
		fscanf(FRV,"%lf",&metal_ISM);
		fscanf(FRV,"%lf",&metal_star_mass);
		fscanf(FRV,"%lf",&metal_star_BL);
		fscanf(FRV,"%lf",&BL);
		fscanf(FRV,"%lf",&tau_V);
		fscanf(FRV,"%lf",&dust);
		fscanf(FRV,"%lf",&SFR);
		fscanf(FRV,"%lf",&Ly);
		fscanf(FRV,"%lf",&SN2);
		fscanf(FRV,"%lf",&SN1a);
		fscanf(FRV,"%lf",&age_star_mass);
		fscanf(FRV,"%lf",&age_star_BL);
		
		/*fluxの読み込み*/
		for(k=1;k<=nbcontinuum; k++)
				fscanf(FRV,"%le",&flux_conti[k]);
		wave_conti[0]=0;
		for(k=1;k<=nblines; k++)
				fscanf(FRV,"%le",&flux_lines[k]);
		wave_lines[0]=0;
		/*波長vs.fluxファイル*/
		if(age==a)
		{
			for(k=1; k<=nbcontinuum; k++)
			{
				for(h=1; h<=nblines; h++)
				{
					if( (wave_conti[k-1] < wave_lines[h] && wave_lines[h]<=wave_conti[k]) && (flux_lines[h] != 0) )
							fprintf(synthspecflux[1],"%lf %le\n",wave_lines[h],flux_lines[h]);
				}
				fprintf(synthspecflux[1],"%lf %le\n",wave_conti[k],flux_conti[k]);
			}
			fprintf(synthspecflux[1],"\n");
			a=a+300;
		}
	}
	
	fclose(FRV);
	fclose(synthspecflux[1]);
	
	return(0);
}
