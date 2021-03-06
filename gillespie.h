#ifndef GILLESPIE_H
#define GILLESPIE_H

#include "structs.h"


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include <gsl/gsl_rng.h>					// random number generator functions
#include <gsl/gsl_randist.h>				// random number distributions
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>


int stochMigration(struct foodweb nicheweb, struct migration, const double[], gsl_rng* rng1, const gsl_rng_type* rng1_T, int, gsl_matrix* Dchoice);
int select_patch(struct migration stochastic,  double, gsl_rng* rng1, int);
int select_patch_random(struct foodweb nicheweb, gsl_rng* rng1);
int select_species(struct foodweb nicheweb, struct migration stochastic, gsl_rng* rng1, int, const double[], int, int);
double choose_time(double, gsl_rng* rng1);
int select_species_random(struct foodweb nicheweb, struct migration stochastic, gsl_rng* rng1);
int stochMigration_version2(struct foodweb nicheweb, struct migration stochastic, const double[], gsl_rng* rng1, const gsl_rng_type* rng1_T, int , gsl_matrix* );
gsl_vector* select_patch_and_species(struct foodweb nicheweb, struct migration stochastic, gsl_vector* , double, int,gsl_rng* rng1, gsl_vector*);

#endif