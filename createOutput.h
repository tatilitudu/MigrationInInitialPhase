#ifndef CREATEOUTPUT_H
#define CREATEOUTPUT_H

#include <gsl/gsl_vector.h>

int createOutputGeneral(struct foodweb nicheweb, struct resource res, struct migration stochastic, char* aims, gsl_vector* robustness, gsl_vector* standardDeviationAll, int L, double mu, double nu, double ymigr, double ymigrDeviation, double migrationEventNumber, double migrationEventNumberDeviation);
int createOutputPatchwise(struct foodweb nicheweb, struct resource res, struct migration stochastic, char* aims, gsl_vector* robustness, gsl_vector* standardDeviationAll, int L, int l);
int createOutputSpeciesNumber(struct foodweb nicheweb, struct resource res, char* aims, double SpeciesNumber[][2], int L, double migrationEventNumber);
int createOutputPatchlink(struct foodweb nicheweb, struct resource res, char* aims, double AllMu[][2], double AllNu[][2], int L, double migrationEventNumber);
int createOutputBiomass(struct foodweb nicheweb, const double y[]);
int createOutputRobustnessPatchwiseEachRun(struct foodweb nicheweb, struct data patchwise[],char* aims, FILE* RobustnessEachRun, int L);

#endif