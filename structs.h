#ifndef STRUCTS_H
#define STRUCTS_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/*
	Die foodweb Struktur enthält die Konsolenparameter S, B, Rnum, Y, T, M, d, x sowie einen Netzwerk Vektor für die Ergebnisse der Nischennetz-Berechnung.
*/
struct foodweb{

	gsl_vector* network;
	gsl_vector* fixpunkte;				// fix0,1,2, fixp0,1,2, testf0,1,2
	gsl_vector* migrPara;
	
	struct simuParams* simuParams;	// Feste Simulationsparameter wie Attackraten etc.
	struct simuMemory* simuMem;		// Fester Simulationsspeicher für Lösen der DGL
	
	
	int S;
	int B;
	int Rnum;

	int Y;
	int T;
	int Tchoice;
	
	double d;
	double x;

	int M;
	int Z;
	
};

struct simuParams{
    const double alpha;		// respiration coefficient
    const double lambda;	// ecologic efficiency
    const double hand;		// handling time
    const double beta;		// intraspecific competition
    const double aij;		// attack rate (same for all species)
};

struct simuMemory{
// muss für jedes neue Nahrungsnetz überschrieben werden, also L mal
    gsl_vector* Masses;		// Massen
    gsl_matrix* A_scaled;	// Adjazenzmatrix mit a und f_ij skaliert
    gsl_vector* tvec;		// Hilfsspeicher einzelnes Habitat
    gsl_vector* rvec;		// Hilfsspeicher einzelnes Habitat 
    gsl_vector* svec;		// Hilfsspeicher einzelnes Habitat
};


struct migration{
  
	gsl_vector* SpeciesNumbers;
	gsl_vector* AllMus;
	gsl_vector* AllNus;
	gsl_vector* Biomass_SpeciesNumbers;
	gsl_vector* Biomass_AllMus;
	gsl_vector* Biomass_AllNus;
	double Bmigr;
	
	gsl_vector *a;
	gsl_vector *c;
	gsl_vector *linkCount;
	
	gsl_vector *aSpecies;
	gsl_vector *cSpecies;
	
};



struct resource{

	double size;
	double growth;

};

struct data{
      gsl_vector* sini;
      gsl_vector* sfini;
      gsl_vector* bini;
      gsl_vector* bfini;
      gsl_vector* robness;
};
      



#endif
