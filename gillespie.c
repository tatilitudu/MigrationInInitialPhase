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

#include "gillespie.h"

#define SEED	123

int stochMigration(struct foodweb nicheweb, struct migration stochastic, const double y[], gsl_rng* rng1, const gsl_rng_type* rng1_T, int migrationEventNumber, gsl_matrix* Dchoice)
{
  int Y = nicheweb.Y;
  int S = nicheweb.S;
  double d = nicheweb.d; 
  int Tchoice = nicheweb.Tchoice;
  int Rnum = nicheweb.Rnum;
  int Z = nicheweb.Z;
  //gsl_vector *network 	= nicheweb.network;		// Inhalt: A+linksA+Y+linksY+Massen+Trophische_Level = (Rnum+S)²+1+Y²+1+(Rnum+S)+S
//   gsl_vector *a	= gsl_vector_calloc(Y);
  double atot;
/*  gsl_vector_view D_view = gsl_vector_subvector(network, (Rnum+S)*(Rnum+S)+1, Y*Y);					// Migrationsmatrix D als Vektor
  gsl_matrix_view ED_mat = gsl_matrix_view_vector(&D_view.vector, Y, Y);								// D als Matrixview
  gsl_matrix *EDmat	 = &ED_mat.matrix;*/	
  
//   gsl_vector *c	= gsl_vector_calloc(Y);
//   gsl_vector *linkCount	= gsl_vector_calloc(Y);
  double ctemp=0;
  double Bmigr = stochastic.Bmigr;
  //printf("Bmigr ist %f\n", Bmigr);
  int l,i;
  double dij 	= pow(10, d);
  
  // Setze c(l) als (Population in Patch l)/(Population in allen Patches) 
  for(l=0;l<Y;l++)
  {
    ctemp = 0;
    for(i=Rnum;i<Rnum+S;i++)
    {
      ctemp += y[l*(Rnum+S)+i];
    }
    gsl_vector_set(stochastic.c,l,ctemp);
    //printf("c[%i] ist %f\n",l,gsl_vector_get(c,l));
  }
  //ctot = gsl_blas_dasum(c);
  gsl_vector_scale(stochastic.c,1/Bmigr);
  
  int linkCountTemp, m;

  
  
  if( Tchoice == 0 )
  {
    for(i=0; i<7;i++)
    {
      gsl_vector_set(nicheweb.migrPara, i, 0);
    }
    printf("Es findet keine Migration statt\n");
  }
  else
  {
    for(l=0;l<Y;l++)
    {
      linkCountTemp=0;
      for(m=0;m<Y;m++)
      {
	if( gsl_matrix_get(Dchoice,l,m)!=0 )
	{
	  linkCountTemp++;
	}
      }
      gsl_vector_set(stochastic.linkCount,l,linkCountTemp);
    }
  
    gsl_vector_memcpy(stochastic.a,stochastic.linkCount);
    gsl_vector_mul(stochastic.a,stochastic.c);
    gsl_vector_scale(stochastic.a,dij);
    //printf("c ist %f\n",gsl_vector_get(c,3));
    //printf("a ist %f\n",gsl_vector_get(a,3));
    atot = gsl_blas_dasum(stochastic.a);
    //printf("atot ist %f\n",atot);
    //printf("Z ist %i\n", Z);
    
    //printf("Berechne Zeitpunkte, zu den migriert werden soll:\n");
    double tau = choose_time(atot, rng1);
    //printf("tau ist %f\n",tau);
    gsl_vector_set(nicheweb.migrPara, 0, tau);

    //printf("\n");
    //r = (double)rand()/INT_MAX;
    //printf("r ist %f\n",r);
  
    //printf("Berechne von welchem Patch aus migriert werden soll\t");
    //printf("r: %f\n",r);
    //printf("atot: %f\n",atot);
    int mu = select_patch(stochastic,  atot, rng1, Y);
    //printf("population ist %f\n",gsl_vector_get(stochastic.Biomass_AllMus, migrationEventNumber));
    //gsl_vector_set(stochastic.AllMus, migrationEventNumber, mu);
    gsl_vector_set(nicheweb.migrPara, 1, mu);
    double biomassAllMus = gsl_vector_get(stochastic.c,mu)*Bmigr;
    double test;
    for(i = 0; i< S ; i++)
    {
      test += y[i+Rnum];
    }
//     printf("über y ist %f\n",test);
//     printf("über c ist %f\n",biomassAllMus);
    
    //gsl_vector_set(stochastic.Biomass_AllMus, migrationEventNumber, biomassAllMus);
    //printf("mu: %i\n",mu);
    int flag=1;
  
    //printf("Berechne in welches Patch migriert werden soll\t\t");
    int nu; 
    while(flag != 0)
    {
      //r1  = (double)rand()/INT_MAX;
      //printf("r1 ist %f\n",r1);
      nu = select_patch_random(nicheweb, rng1);
      if(nu!= mu  && gsl_matrix_get(Dchoice, nu, mu)!=0)
      {
	gsl_vector_set(nicheweb.migrPara, 2, nu);
	//gsl_vector_set(stochastic.AllNus, migrationEventNumber, nu); 
	flag = 0;
      }
    }
    //printf("nu: %i\n", nu);
    
    
    int SpeciesNumber;
    //r2 = (double)rand()/INT_MAX;
    //printf("r2 ist %f\n",r2);
    //int Choice = 0;
    SpeciesNumber = select_species(nicheweb, stochastic, rng1, 0, y, migrationEventNumber, mu);
    //gsl_vector_set(stochastic.SpeciesNumbers, migrationEventNumber, SpeciesNumber);
    //gsl_vector_set(stochastic.Biomass_SpeciesNumbers, migrationEventNumber, y[(S+Rnum)*mu+(SpeciesNumber+Rnum)]);
    gsl_vector_set(nicheweb.migrPara, 3, SpeciesNumber);
    
    //printf("SpeciesNumber: %i\n\n", SpeciesNumber);
    //printf("Population dieser Spezies ist %f\n",y[SpeciesNumber+Rnum]);
    
    if(SpeciesNumber>S)
    {
      printf("\n\nFehler!!! SpeciesNumber>S \n\n");
      
    }
    
//     gsl_vector_free(a);
  }
  

  
//   gsl_vector_free(c);
//   gsl_vector_free(linkCount);
  
  
  
  return 0;
}

double choose_time(double atot, gsl_rng* rng1)
{
  double tau = 0;
  double r = gsl_rng_uniform_pos(rng1);
  //printf("Zufallszahl ist %f\n",r);
  //printf("atot ist %f\n",atot);
  if( atot>0 )
  {
    
    // rand liefert zufällige Zahl zwischen 0 und INT_MAX
    tau = -log(r)/ atot;
  }
    
  return tau;
}


int select_patch(struct migration stochastic,  double atot, gsl_rng* rng1, int Y)
{
  double r = gsl_rng_uniform_pos(rng1);
  int i;
  int mu;
  double sum=0;
  
  r = r*atot;
  //printf("r*atot in select_patch ist %f\n",r);
  for(i=0;i<Y;i++)
  {
//     printf("a[%i] = %f",i,gsl_vector_get(stochastic.a,i));
    sum += gsl_vector_get(stochastic.a,i);
    if( r < sum )
    {
      mu = i;
      break;
    }
  }
//   if(whichPatch == 1)
//   {
//     gsl_vector_set(stochastic.Biomass_AllNus, migrationEventNumber, gsl_vector_get(a,mu));
//   }
//   else if(whichPatch == 0)
//   {
//     
//     gsl_vector_set(stochastic.Biomass_AllMus, migrationEventNumber, gsl_vector_get(a,mu));
//     //printf("population1 ist %f\n",gsl_vector_get(stochastic.Biomass_AllMus, migrationEventNumber));
//   }
    
  return mu;
}


int select_patch_random(struct foodweb nicheweb, gsl_rng* rng1)
{
  double r = gsl_rng_uniform_pos(rng1);
  
  int Y = nicheweb.Y;
  int patchNumber;
  float patchNumberFloat;
  //r = 0;
  patchNumberFloat = r*(Y);
//   printf("patchNumberFloat ist %f\n", patchNumberFloat);
  
  patchNumber = (int)patchNumberFloat;
  if(patchNumber>3)
  {
    patchNumber = 3;
//     printf("nu ist %i\n",patchNumber);
//     printf("patchNumberFloat ist %f\n",patchNumberFloat);
//     printf("Zufallszahl ist %f\n", r);
  }
//   if(patchNumberFloat>0) patchNumber = (int)(patchNumberFloat + 0.5);
// 	
//   else patchNumber =  (int)(patchNumberFloat - 0.5);
  
//   printf("patchNumber ist %i\n", patchNumber);
  
  return patchNumber;
}


int select_species(struct foodweb nicheweb, struct migration stochastic, gsl_rng* rng1, int Choice, const double y[], int migrationEventNumber, int mu)
{
  int S = nicheweb.S;
  int Y = nicheweb.Y;
  int Rnum = nicheweb.Rnum;
//   gsl_vector *network = nicheweb.network;
//   gsl_vector_view M_vec = gsl_vector_subvector(network, ((Rnum+S)*(Rnum+S))+1+(Y*Y)+1, (Rnum+S));	// Massenvektor
//   gsl_vector *Mvec = &M_vec.vector;
  
  
  
  int i;
  int SpeciesNumber;
  double sum = 0;
  double atot;
  double r = gsl_rng_uniform_pos(rng1);
//   gsl_vector *a = gsl_vector_calloc(S);
//   gsl_vector *c = gsl_vector_calloc(S);
  
  
  //printf("Berechne, welche Spezies migrieren darf\t\t\t");
  
  for(i = 0; i< S ;i++)
  {
    gsl_vector_set(stochastic.aSpecies,i,y[mu*(Rnum+S)+(Rnum+i)]);
  }
  //printf("Eintrag 5 von y ist %f\n", y[5]);
  // Nur Abhängigkeit, wie groß Population ist == 0; zusätzlich Massen == 1
  if(Choice == 0)
  {
    atot = gsl_blas_dasum(stochastic.aSpecies);
    gsl_vector_scale(stochastic.aSpecies, 1/atot);
  }
  else if(Choice == 1)
  {
    double ctot;
    for(i = 0; i < S; i++ )
    {
      gsl_vector_set(stochastic.cSpecies, i , gsl_vector_get(nicheweb.simuMem->Masses,i+Rnum));
    }
    ctot = gsl_blas_dasum(stochastic.cSpecies);
    gsl_vector_scale(stochastic.cSpecies,1/ctot);
    gsl_vector_mul(stochastic.aSpecies,stochastic.cSpecies);
    atot = gsl_blas_dasum(stochastic.aSpecies);
  }
    
  // Abhängigkeit
  
  //printf("atot ist %f\n",atot);
  
  //r = r*atot;

  //printf("r*atot ist %f\n",r);
  
  for(i=0;i<S;i++)
  {
    sum += gsl_vector_get(stochastic.aSpecies,i);
    if( r < sum )
    {
      SpeciesNumber = i;
      break;
    }
  }
//   double biomassSpecies = 1/atot*gsl_vector_get(a,SpeciesNumber);
//   gsl_vector_set(stochastic.Biomass_SpeciesNumbers, migrationEventNumber, biomassSpecies);

//   gsl_vector_free(stochastic.aSpecies);
//   gsl_vector_free(c);
  
  return SpeciesNumber;
  
}

int stochMigration_version2(struct foodweb nicheweb, struct migration stochastic, const double y[], gsl_rng* rng1, const gsl_rng_type* rng1_T, int migrationEventNumber, gsl_matrix* Dchoice)
{
  double d = nicheweb.d; 
  int Tchoice = nicheweb.Tchoice;
  int Rnum = nicheweb.Rnum;
  int Z = nicheweb.Z;
  int S = nicheweb.S;
  int Y = nicheweb.Y;
  //gsl_vector *network 	= nicheweb.network;		// Inhalt: A+linksA+Y+linksY+Massen+Trophische_Level = (Rnum+S)²+1+Y²+1+(Rnum+S)+S
  
/*  gsl_vector_view D_view = gsl_vector_subvector(network, (Rnum+S)*(Rnum+S)+1, Y*Y);					// Migrationsmatrix D als Vektor
  gsl_matrix_view ED_mat = gsl_matrix_view_vector(&D_view.vector, Y, Y);								// D als Matrixview
  gsl_matrix *EDmat	 = &ED_mat.matrix;*/	
  gsl_vector *selection = gsl_vector_calloc(2);
  gsl_vector *c	= gsl_vector_calloc(S*Y);
  gsl_vector *linkCount	= gsl_vector_calloc(Y);
  //double ctemp=0;
  double Bmigr = stochastic.Bmigr;
  int l,i;
  double dij 	= pow(10, d);
  //printf("\n");
  gsl_vector *a	= gsl_vector_calloc(Y*(S));
    double atot;
  // Setze c(l) als (Population in Patch l)/(Population in allen Patches) 
  for(l=0;l<Y;l++)
  {
    for(i=Rnum;i<Rnum+S;i++)
    {
      gsl_vector_set(c,l*S+(i-Rnum),y[l*(Rnum+S)+i]);
//       printf("Index ist %i und %i\n",l,i);
//       printf("y ist %f\t", y[l*(Rnum+S)+i]);
//       printf("c ist %f\t", gsl_vector_get(c,l*S+(i-Rnum)));
    }
//     printf("\n\n");
  }
//   printf("Raus aus der Schleife\n");
  gsl_vector_scale(c,1/Bmigr);
  
//   for(l = 0; l<Y;l++)
//   {
//     for(i=Rnum;i<Rnum+S;i++)
//     {
//     //if(gsl_vector_get(c,l)!=0)
//     //{
//       printf("c ist %f\n",gsl_vector_get(c,l));
//     }
//     //}
//   }

  int linkCountTemp, m;
  
 
  if( Tchoice == 0 )
  {
    for(i=0; i<3*Z+3;i++)
    {
      gsl_vector_set(nicheweb.migrPara, i, 0);
    }
    printf("Es findet keine Migration statt\n");
  }
  else
  {
    for(l=0;l<Y;l++)
    {
      linkCountTemp=0;
      for(m=0;m<Y;m++)
      {
	if( gsl_matrix_get(Dchoice,l,m)!=0 )
	{
	  linkCountTemp++;
	}
      }
      for(i = 0; i< S; i++)
      {
	gsl_vector_set(a, l*(S)+i, linkCountTemp);
	//printf("Ergebnis ist %f\n", gsl_vector_get(a,l*(S)+i));
      }
      //gsl_vector_set(linkCount,l,linkCountTemp);
    }
    //gsl_vector_scale(linkCount,1/linkCountTot);
  
  
    
    //double r,r1,r2;
  
    
    //gsl_vector_memcpy(a,linkCount);
    gsl_vector_mul(a,c);
    gsl_vector_scale(a,dij);
  
    atot = gsl_blas_dasum(a);
    

    double tau = choose_time(atot, rng1);
    gsl_vector_set(nicheweb.migrPara, 0, tau);

    //printf("tau ist %f\n",tau);
    
    selection = select_patch_and_species(nicheweb, stochastic, a, atot, migrationEventNumber, rng1, selection);
    int mu = gsl_vector_get(selection,0);  
    int SpeciesNumber = gsl_vector_get(selection,1);
//     printf("mu ist draußen %i\n",mu);
//     printf("Spezies ist draußen %i\n",SpeciesNumber);
    
    
    //printf("population ist %f\n",gsl_vector_get(stochastic.Biomass_AllMus, migrationEventNumber));
    gsl_vector_set(stochastic.AllMus, migrationEventNumber, mu);
    double biomassAllMus = 0;
    for(i = 0; i< S; i++)
    {
      biomassAllMus += y[(S+Rnum)*mu+(i+Rnum)];
    }
    
    gsl_vector_set(stochastic.Biomass_AllMus, migrationEventNumber, biomassAllMus);
    //printf("y von Spezies %i auf patch %i ist %f\n",SpeciesNumber,mu,y[(S+Rnum)*mu+(SpeciesNumber+Rnum)]);
    gsl_vector_set(stochastic.Biomass_SpeciesNumbers, migrationEventNumber, y[(S+Rnum)*mu+(SpeciesNumber+Rnum)]);
    gsl_vector_set(nicheweb.migrPara, 1, mu);
    //printf("mu: %i\n",mu);
    
    gsl_vector_set(stochastic.SpeciesNumbers, migrationEventNumber, SpeciesNumber);
    gsl_vector_set(nicheweb.migrPara, 3, SpeciesNumber);
    
    
    if(SpeciesNumber>S)
    {
      printf("\n\nFehler!!! SpeciesNumber>S \n\n");
      
    }
    
    int flag=1;
  
    //printf("Berechne in welches Patch migriert werden soll\t\t");
    int nu; 
    while(flag != 0)
    {
      //r1  = (double)rand()/INT_MAX;
      //printf("r1 ist %f\n",r1);
      nu = select_patch_random(nicheweb, rng1);
      
      if(nu!= mu  && gsl_matrix_get(Dchoice, nu, mu)!=0)
      {
	gsl_vector_set(nicheweb.migrPara, 2, nu);
	gsl_vector_set(stochastic.AllNus, migrationEventNumber, nu); 
	flag = 0;
      }
    }
    //printf("nu: %i\n", nu);
    
   
    
    
    gsl_vector_free(a);
  }
  

  
  gsl_vector_free(c);
  gsl_vector_free(linkCount);
  gsl_vector_free(selection);
  
  
  
  return 0;
}


gsl_vector* select_patch_and_species(struct foodweb nicheweb, struct migration stochastic, gsl_vector* a, double atot, int migrationEventNumber,gsl_rng* rng1, gsl_vector* selection)
{
  //int Rnum = nicheweb.Rnum;
  int S = nicheweb.S;
  int Y = nicheweb.Y;
  double r = gsl_rng_uniform_pos(rng1);
  int i,j;
  double sum = 0;
  
  bool stop = false;
  r = r*atot;
//   printf("r ist %f\n",r);
  
  for(i = 0 ; i< Y && !stop; i++)
  {
    for(j = 0; j < S && !stop; j++)
    {
      sum += gsl_vector_get(a,(S)*i+j);
//       printf("a ist %f\t",gsl_vector_get(a,(S)*i+j));
//       printf("sum ist %f\t",sum);
      if(r< sum)
      {
// 	printf("\n mu ist %i\n",i);
// 	printf("\n Spezies ist %i\n",j);
	gsl_vector_set(selection, 0, i);
	gsl_vector_set(selection, 1, j);
	stop = true;
	
	break;
      }
    }
  }
//   printf("raus\n");
  return selection;
}

int select_species_random(struct foodweb nicheweb, struct migration stochastic, gsl_rng* rng1)
{
  double r = gsl_rng_uniform_pos(rng1);
  
  int S = nicheweb.S;
  int SpeciesNumber;
  float SpeciesNumberFloat;
  
  
  //printf("Berechne, welche Spezies migrieren darf\t\t\t");
  //printf("Zufallszahl ist %f\n",r);
  
  SpeciesNumberFloat = r*(S);
//   printf("SpeciesNumberFloat ist %f\n", SpeciesNumberFloat);
  
  SpeciesNumber = (int)SpeciesNumberFloat;
//   if(SpeciesNumberFloat>0) SpeciesNumber = (int)(SpeciesNumberFloat + 0.5);
// 	
//   else SpeciesNumber =  (int)(SpeciesNumberFloat - 0.5);
  
//   printf("SpeciesNumber ist %i\n", SpeciesNumber);
  
  
  
  return SpeciesNumber;
}