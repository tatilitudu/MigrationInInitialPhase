/* 2015-07-21 14:50:40 

	Projekt: Nischennetz mit Migration mit Robustness Analyse

	Quelldatei: main.c zu Nischennetz_mit_Migration_1.0
	Programmiersprache: C
	Autoren: S. J. Plitzko, M.Hamm
	Version 1.0

###############################################################################################################################
Beschreibung:
Dieses Programm simuliert ein über mehrere Patches ausgedehntes Nahrungsnetz mit dem Nischenmodell. Die Populationsdynamik wird mit einer Holling Typ II Funktion modelliert. Beim Aufrufen können über die Konsole folgende Parameter bestimmt werden:

Anzahl Spezies S
Anzahl basale Spezies B
Allometrie Koeffizient x
Migrationsstärke d
Topologie des Lebensrausm (siehe Set_Topology())
Anzahl der statistischen Wiederholungen L
Anzahl der Patches Y
Migration ja oder nein "M"

Kompilieren mit: gcc -o V1_15_07_13 main.c getargs.c topology.c mignicheweb.c evolveweb.c holling2.c robustness.c -lm -lgsl -lgslcblas -Wall

Eingabe in der Konsole: ./NAME -S X -B X -T X -L X -Y X -d X.X -x X.X -M "X" -R X.X	mit X bzw. X.X numerischer Wert der Variablen (Reihenfolge egal)

Das Ergebnis der Simulation wird in der Datei xxx.out gespeichert. Dabei sind die Daten wie folgt angeordnet (alles eine Zeile): 
S		B		M		x		Y		dpow	T
Rob		Perlok	Perges
Si_ges	Si_TL1	Si_TL2	Si_TL3	Si_TL4	Si_TL>4		Spezies initial
Sf_ges	Sf_TL1	Sf_TL2	Sf_TL3	Sf_TL4	Sf_TL>4		Spezies final
Bi_ges	Bi_TL1	Bi_TL2	Bi_TL3	Bi_TL4	Bi_TL>4		Biomass initial
Bf_ges	Bf_TL1	Bf_TL2	Bf_TL3	Bf_TL4	Bf_TL>4		Biomass final
Sh_ges	Sh_TL1	Sh_TL2	Sh_TL3	Sh_TL4	Sh_TL>4		Spezies Hub
Bh_ges	Bh_TL1	Bh_TL2	Bh_TL3	Bh_TL4	Bh_TL>4		Biomass Hub
Ss_ges	Ss_TL1	Ss_TL2	Ss_TL3	Ss_TL4	Ss_TL>4		
Bs_ges	Bs_TL1	Bs_TL2	Bs_TL3	Bs_TL4	Bs_TL>4
1mit2	2mit3	3mit1
Fixp0	Fixp1	Fixp2	Fixp3	Fixp4	Fixp5	Fixp6	Fixp7

###############################################################################################################################
*/

//--Standardbibliotheken--------------------------------------------------------------

#include <string.h>						// string modification functions
#include <time.h>						// time functions
#include <math.h>						// math functions
#include <stdio.h>						// output functions
#include <stdlib.h>						// standard
#include <limits.h>
#include <gsl/gsl_rng.h>				// random number generator functions
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>			// random number distributions
#include <gsl/gsl_sort_vector.h>		// vector sorting functions
#include <gsl/gsl_odeiv.h>				// differential equation solver
#include <gsl/gsl_errno.h>				// errorhandler


//--Selbstdefinierte Funktionen------------------------------------------------------
#include "structs.h"					// foodweb Struktur

#include "getargs.h"					// getArgs
#include "topology.h"					// SetTopology
#include "mignicheweb.h"				// SetNicheNetwork, LinkElements, CountLinks
#include "evolveweb_neu.h"					// DGL lösen
#include "holling2_festeMigrmenge_neu.h"					// DGL lösen			
#include "robustness.h"					// Robustness Analyse
#include "StandardDeviation.h"				// Standardabweichung berechnen
#include "createOutput.h"				// Output erzeugen

//--Verzeichnis für Ergebnisse-------------------------------------------------------
#define ORT "/home/tatjana/Arbeitsfläche/MichaelasProgramm/stochastischeMigration/Migration_in_Anfangsphase/Erste_Versuche/"
#define ORT2 "/home/tatjana/Arbeitsfläche/MichaelasProgramm/stochastischeMigration/Migration_in_Anfangsphase/Erste_Versuche/"
//++START++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int main(int argc, char** argv)
{	

//--Foodweb Struktur mit Standardwerten aufstellen------------------------------------------------------------------------------------------
	struct simuParams simParams 	= {0.3, 0.65, 0.35, 0.5, 6.0};			// Diese Parameter sind konstant
	struct simuMemory simMem 	= {NULL, NULL, NULL, NULL, NULL};		// Größe der Vektoren liegt noch nicht fest	

	gsl_vector* fixpunkte	= gsl_vector_calloc(9);
 	

	struct foodweb nicheweb	= {NULL, fixpunkte, NULL,&simParams, &simMem, 18, 3, 1, 5, 0, 0, -7., 0.0, 0, 1};		// Reihenfolge: network, fxpkt, migrPara, AllMus, AllNus, S, B, Rnum, Y, T, Tchoice, d, x, M, Z
	
	struct migration stochastic = {NULL, NULL, NULL, NULL, NULL, NULL, 0.00001, NULL, NULL, NULL, NULL, NULL};
	
	struct resource res = {500.0, 0.0};											// Resource: Größe, Wachstum
	
	
	
	
//--Konsoleneingabe-------------------------------------------------------------------------------------------------------------------------
	
	int L = 5;	// Statistik
	int i = 0,j;	// Counter
	char aims6[255] = ORT;
	FILE* RobustnessEachRun;
	
	int checksum = getArgs(argc, argv, &(nicheweb.S), &(nicheweb.B), &(nicheweb.T), &(nicheweb.d), &L, &(nicheweb.Y), &(nicheweb.x), &(nicheweb.M), &(res.size), &(nicheweb.Z), &(stochastic.Bmigr));	

		if (checksum != 11 && checksum!=(int)(argc-1)/2) 	// Alles gesetzt?									
		 {	
			printf("Bitte gültige Eingabe für Parameter machen!\nProgramm wird beendet.\n");		
			return(0);		
		 }

/*	int length			= ((nicheweb.Rnum+nicheweb.S)*(nicheweb.S+nicheweb.Rnum)+1+nicheweb.Y*nicheweb.Y+1+(nicheweb.Rnum+nicheweb.S)+nicheweb.S+1);	// Länge des Rückabewerts
	nicheweb.network = gsl_vector_calloc(length);	 
*/	
	// Speicher belegen, nachdem die Größe des Systems durch die Konsoleneingabe bekannt ist	
	CallocFoodwebMem(&nicheweb); 
	CallocStochasticMem(&stochastic, nicheweb.Y, nicheweb.S);
	
	printf("Z = %i\n",nicheweb.Z);
	nicheweb.migrPara = gsl_vector_calloc(7); // Reihenfolge: tau, mu, nu, SpeciesNumber, momentanes t, ymigr, migrationEventNumber	 
// 	stochastic.SpeciesNumbers = gsl_vector_calloc(nicheweb.Z);
// 	stochastic.AllMus = gsl_vector_calloc(nicheweb.Z);
// 	stochastic.AllNus = gsl_vector_calloc(nicheweb.Z);
// 	stochastic.Biomass_SpeciesNumbers = gsl_vector_calloc(nicheweb.Z);
// 	stochastic.Biomass_AllMus = gsl_vector_calloc(nicheweb.Z);
// 	stochastic.Biomass_AllNus = gsl_vector_calloc(nicheweb.Z);
	
//--Zufallszahlengenerator initialisieren--------------------------------------------------------------------------------

		const gsl_rng_type *rng1_T;											// ****
		gsl_rng *rng1;   													// initialize random number generator
		gsl_rng_env_setup();   												// ermöglicht Konsolenparameter
		rng1_T = gsl_rng_default;   										// default random number generator (so called mt19937)
		gsl_rng_default_seed = 0;											// default seed for rng
// 		gsl_rng_default_seed = ((unsigned)time(NULL));						// random starting seed for rng
		rng1 = gsl_rng_alloc(rng1_T);
		
		
		
		
//--Struct initialisieren für patchweise Ausgabe----------------------------------------------------------------------------------------		 

	struct data patchwise[nicheweb.Y];
	for(i=0; i<nicheweb.Y; i++)
	{
	  gsl_vector* sini  = gsl_vector_calloc(6);
	  gsl_vector* sfini  = gsl_vector_calloc(6);
	  gsl_vector* bini  = gsl_vector_calloc(6);
	  gsl_vector* bfini  = gsl_vector_calloc(6);
	  gsl_vector* robness  = gsl_vector_calloc(2);
	  
	  //struct data tempo = {sini,sfini,bini,bfini,robness};
	  struct data temp = {sini,sfini,bini,bfini,robness};
	  patchwise[i] = temp;
	}
	
	//printf("test");
//--Initialisierungen---------------------------------------------------------------------------------------------------
	nicheweb.Tchoice = nicheweb.T;
	nicheweb.T = 0;
	nicheweb.d = nicheweb.d/10;
	printf("d ist %f\n",nicheweb.d);
	res.size = res.size/10;
	stochastic.Bmigr = 2.5*stochastic.Bmigr*pow(10,nicheweb.d);
	nicheweb.Z = pow(10,nicheweb.Z);
	printf("Bmigr ist %f\n",stochastic.Bmigr);
	printf("Z ist %i\n",nicheweb.Z);
	printf("x ist %f\n",nicheweb.x);
	//int len	= ((nicheweb.Rnum+nicheweb.S)*(nicheweb.S+nicheweb.Rnum)+1+nicheweb.Y*nicheweb.Y+1+(nicheweb.Rnum+nicheweb.S)+nicheweb.S);	// Länge des Rückabewerts

	gsl_vector *populationFIN 	= gsl_vector_calloc((nicheweb.Rnum + nicheweb.S)*(nicheweb.Y)*5 + (nicheweb.S) + 3);				// Gleiche Länge wie Rückgabe von evolveNetwork
	gsl_vector *robustness		= gsl_vector_calloc(63);
	gsl_vector *resultEvolveWeb	= gsl_vector_calloc((nicheweb.Rnum+nicheweb.S)*nicheweb.Y*5 + 3 + nicheweb.S); 				// y[Simulation], y0, ymax, ymin, yavg, fixp, TL
	gsl_vector *resultRobustness 	= gsl_vector_calloc(63);
	gsl_matrix *D 			= gsl_matrix_calloc(nicheweb.Y,nicheweb.Y);
	gsl_matrix* Dchoice		= gsl_matrix_calloc(nicheweb.Y,nicheweb.Y);
	
	gsl_vector *robustnesstemp	= gsl_vector_calloc(63);
	gsl_vector *meanSquOfDataAll 	= gsl_vector_calloc(63);
	gsl_vector *meanSquOfDataAlltemp = gsl_vector_calloc(63);
	
	gsl_vector *standardDeviationAll = gsl_vector_calloc(63);
	
	gsl_vector *meanOfData	= gsl_vector_calloc((6*4+2)*nicheweb.Y);
	gsl_vector *meanOfDatatemp = gsl_vector_calloc((6*4+2)*nicheweb.Y);
	gsl_vector *meanSquOfData = gsl_vector_calloc((6*4+2)*nicheweb.Y);
	
	gsl_vector *standardDeviation = gsl_vector_calloc((6*4+2)*nicheweb.Y);
	
	gsl_vector_set_zero(robustness);
	gsl_vector_set_zero(meanOfData);
	gsl_vector_set_zero(meanSquOfData); 
	gsl_vector_set_zero(nicheweb.migrPara);
	gsl_vector_set_zero(meanSquOfDataAll);
	
// 	double SpeciesNumber[L*nicheweb.Z][2]; 
// 	double AllMu[L*nicheweb.Z][2];
// 	double AllNu[L*nicheweb.Z][2];
	
	double ymigr = 0;
	double mu = 0;
	double nu = 0;
	double ymigrtemp;
	double ymigrSqu = 0;
	double ymigrDeviation;
	double migrationEventNumber = 0;
	double migrationEventNumbertemp;
	double migrationEventNumberSqu = 0.0;
	double migrationEventNumberDeviation;
	int lastMigrationEventNumber = 0;

//--Simulation---------------------------------------------------------------------------------------------------------------------
	//SetTopology(nicheweb.Y, nicheweb.T, D);	
	SetTopology(nicheweb.Y, nicheweb.Tchoice, Dchoice);	
// 	for(i = 0; i<nicheweb.Y; i++)
// 	{
// 	  for(j = 0 ; j<nicheweb.Y; j++)
// 	  {
// 	    printf("%f\t",gsl_matrix_get(Dchoice,i,j));
// 	  }
// 	  printf("\n");
// 	}

	for(i = 0; i < L; i++)																							
	 { 	
// 		const gsl_rng_type *rng1_T;											// ****
// 		gsl_rng *rng1;   													// initialize random number generator
// 		gsl_rng_env_setup();   												// ermöglicht Konsolenparameter
// 		rng1_T = gsl_rng_default;   										// default random number generator (so called mt19937)
// 		gsl_rng_default_seed = 0;											// default seed for rng
// 		//gsl_rng_default_seed = ((unsigned)time(NULL));						// random starting seed for rng
// 		rng1 = gsl_rng_alloc(rng1_T);
		
		
		printf("\nStarte Durchlauf L = %i\n", i);
//--Starte Simulation-----------------------------------------------------------------------------------------------			
		SetNicheNetwork(nicheweb, stochastic, res, rng1, rng1_T, D);
		gsl_vector_set_zero(resultEvolveWeb);
		populationFIN	 = EvolveNetwork(nicheweb, stochastic, rng1, rng1_T, Dchoice, resultEvolveWeb);
			
		gsl_vector_set_zero(resultRobustness);										
		gsl_vector_memcpy(robustnesstemp, EvaluateRobustness(populationFIN, nicheweb, patchwise, resultRobustness));	// Robustness Analyse
		
		
//--Standardabweichung für Mittelung vorbereiten-----------------------------------------------------------------------------------------		
		determineMean(robustnesstemp, 63, robustness);
		determineMeanSqu(robustnesstemp, 63, meanSquOfDataAll);
		
//--Ausgabewerte----------------------------------------------------------------------------------------------------------		
		ymigrtemp = gsl_vector_get(nicheweb.migrPara, 5);
		migrationEventNumbertemp = gsl_vector_get(nicheweb.migrPara, 6);
// 		for(int j= 0; j<migrationEventNumbertemp; j++)
// 		{
// 		  AllMu[lastMigrationEventNumber+j][0] = gsl_vector_get(stochastic.AllMus, j);
// 		  AllNu[lastMigrationEventNumber+j][0] = gsl_vector_get(stochastic.AllNus, j);
// 		  SpeciesNumber[lastMigrationEventNumber+j][0] = gsl_vector_get(stochastic.SpeciesNumbers,j);
// 		  
// 		  AllMu[lastMigrationEventNumber+j][1] = gsl_vector_get(stochastic.Biomass_AllMus, j);
// 		  AllNu[lastMigrationEventNumber+j][1] = gsl_vector_get(stochastic.Biomass_AllNus, j);
// 		  SpeciesNumber[lastMigrationEventNumber+j][1] = gsl_vector_get(stochastic.Biomass_SpeciesNumbers,j);
// 		}
		lastMigrationEventNumber += migrationEventNumbertemp;
		
		//printf("SpeciesNumber ist %f\n",SpeciesNumber[i]);
		ymigr += ymigrtemp;
		migrationEventNumber += migrationEventNumbertemp;
		
		ymigrSqu += (ymigrtemp*ymigrtemp); 
		migrationEventNumberSqu = (migrationEventNumberSqu*(i)+(migrationEventNumbertemp*migrationEventNumbertemp))/(i+1); 
// 		printf("Additionsteil ist %f\n",(migrationEventNumberSqu*(i)+(migrationEventNumbertemp*migrationEventNumbertemp))/(i+1));
//--Mittelwert und Vorbereitungen für Standardabweichung für die patchweise Ausgabe berechnen--------------------------------		
		linkElements(patchwise, nicheweb.Y, meanOfDatatemp);
		
		determineMean(meanOfDatatemp, (6*4+2)*nicheweb.Y, meanOfData);
		determineMeanSqu(meanOfDatatemp, (6*4+2)*nicheweb.Y, meanSquOfData);
		
		createOutputRobustnessPatchwiseEachRun(nicheweb,patchwise,aims6,RobustnessEachRun,i);
		
		printf("\nBeende Durchlauf L = %i\n", i);
	 }


//-- Standardabweichung berechnen--------------------------------------------------------------------------------------

	ymigrSqu = ymigrSqu/L;
	ymigr = ymigr/L;
	ymigrDeviation = sqrt(ymigrSqu - ymigr*ymigr);
	
	//migrationEventNumberSqu = migrationEventNumberSqu/L;
	migrationEventNumber = migrationEventNumber/L;
	migrationEventNumberDeviation = sqrt(migrationEventNumberSqu - migrationEventNumber*migrationEventNumber);
// 	printf("migrationEventNumber ist %f\n",migrationEventNumber);
// 	printf("migrationEventNumberSqu ist %f\n",migrationEventNumberSqu);
// 	printf("migrationEventNumberDeviation ist %f\n",migrationEventNumberDeviation);
// 	
//-- Für patchweise Ausgabe-------------------------------------------------------------------------------------------
	standardDeviation = determineStandardDeviation((6*4+2)*nicheweb.Y, meanOfData, meanSquOfData, L, standardDeviation);
	standardDeviationAll = determineStandardDeviation(63, robustness, meanSquOfDataAll, L, standardDeviationAll);
	 //printf("der 3. Eintrag in standardDeviationAll ist %f\n", gsl_vector_get(standardDeviationAll,3));
// 	 printf("S ist %f\n", gsl_vector_get(robustness,3));
// 	 printf("Standardabweichung von S ist %f\n", gsl_vector_get(standardDeviationAll,3));
// 	 printf("meanOfDataSqu ist %f\n", gsl_vector_get(meanOfDataSquAll,3));
// 	 printf("meanSquOfData ist %f\n", gsl_vector_get(meanSquOfDataAll,3));
	 
	 
	 
	printf("L=%i\tspeciesini=%f\tspeciesfinal=%f\n", L, gsl_vector_get(robustness, 3)/L, gsl_vector_get(robustness, 9)/L);
	

	
	
//--Abspeichern in File-------------------------------------------------------------------------------------	
	char aims[255] = ORT;

	
	createOutputGeneral(nicheweb, res, stochastic, aims, robustness, standardDeviationAll, L, mu, nu, ymigr, ymigrDeviation, migrationEventNumber, migrationEventNumberDeviation);													// Datei schließen

	
    
	
//--Daten patchweise abspeichern----------------------------------------------------------------------	
// 	printf("population ist %f\n",gsl_vector_get(stochastic.Biomass_AllMus,0));
     
     for(int l = 0 ; l< nicheweb.Y; l++)
     {
       //char name[100];
       char aims2[255] = ORT2;
       
       createOutputPatchwise(nicheweb, res, stochastic, aims2, meanOfData, standardDeviation, L, l);
       
      }
      
//       if(nicheweb.Tchoice != 0)
//       {
//--Ausgewählte Spezies rausschreiben, die migrieren darf---------------------------------------------------------------------------
	
// 	char aims3[255] = ORT;
// 	
// 	//createOutputSpeciesNumber(nicheweb, res, aims3, SpeciesNumber, L, migrationEventNumber);
// 
//       
// //--Ausgewählte Verbindung rausschreiben, über die migriert werden darf---------------------------------------------------------------------------
// 	
// 	char aims4[255] = ORT;
// 	
// // 	createOutputPatchlink(nicheweb, res, aims4, AllMu, AllNu, L, migrationEventNumber);
//       }
      
      printf("\nSimulation abgespeichert\n\n");
	
//--free----------------------------------------------------------------------------------------------------------------  
	free(nicheweb.network);
	
	gsl_vector_free(fixpunkte);
	
	for(i=0; i<nicheweb.Y; i++)
	{
	  gsl_vector_free(patchwise[i].sini);
	  gsl_vector_free(patchwise[i].sfini);
	  gsl_vector_free(patchwise[i].bini);
	  gsl_vector_free(patchwise[i].bfini);
	  gsl_vector_free(patchwise[i].robness);

	}
	
	FreeFoodwebMem(&nicheweb);	// eigene Funktion
	FreeStochasticMem(&stochastic);
	gsl_vector_free(nicheweb.migrPara);
// 	gsl_vector_free(stochastic.AllMus);
// 	gsl_vector_free(stochastic.AllNus);
// 	gsl_vector_free(stochastic.SpeciesNumbers);
// 	gsl_vector_free(stochastic.Biomass_AllMus);
// 	gsl_vector_free(stochastic.Biomass_AllNus);
// 	gsl_vector_free(stochastic.Biomass_SpeciesNumbers);
	gsl_vector_free(populationFIN);
	gsl_vector_free(robustness);	
	gsl_vector_free(meanOfData);
	gsl_vector_free(meanOfDatatemp);
	gsl_vector_free(meanSquOfData);
	gsl_vector_free(standardDeviation);
	gsl_vector_free(standardDeviationAll);
	gsl_vector_free(meanSquOfDataAll);
	gsl_vector_free(meanSquOfDataAlltemp);
	gsl_matrix_free(D);
	gsl_matrix_free(Dchoice);
	//gsl_vector_free(resultEvolveWeb);
	gsl_vector_free(robustnesstemp);
	gsl_rng_free(rng1);
	
	return(0);

}




















