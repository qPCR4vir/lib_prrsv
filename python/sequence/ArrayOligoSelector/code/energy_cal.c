#ifndef ENERGY_CAL
#define ENERGY_CAL

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "data.h"

#define MATCH 0 //middle line in blastoutput
#define PLUS 1 //first line in blastoutput
#define MINUS 2  //third line in blastoutput

#define MAX_LOOP 30
#define R 1.9875*0.001 // all energy parameters are in kcal/mol unit
#define bool int

#endif

double WC_Matrix[5][5][5][5];
double loop[MAX_LOOP+1];
double bulge[MAX_LOOP+1];
double F_function[6];  

void Initail_Parameter();

double Perfect_Energy (int start, int end, int** DATA)
{
  int i;
  double energy=0;

  if (start>= end)
    {
      return 0;
    }
  //printf ("%d %d \n", start, end);

  for (i=start; i<end; i++)
    {  
      //debugging
      //printf("%f\n",  WC_Matrix[DATA[PLUS][i]][DATA[PLUS][i+1]][DATA[MINUS][i]][DATA[MINUS][i+1]] );
      energy += WC_Matrix[DATA[PLUS][i]][DATA[PLUS][i+1]][DATA[MINUS][i]][DATA[MINUS][i+1]];
    }
  //debugging
  //  printf("wc %f\n", energy);
  return energy;

}

int Find_nogap (int * data, int start, int end)
{
  int i;
  
  if (start< end)
    for (i=start; i<end;i++)
      {
	if (data[i] != 0)
	  return i;
      }
  else if (start > end)
    for (i=start; i>end;i--)
      {
	if (data[i] != 0)
	  return i;
      }
  return -1;
}

/*http://bioinfo.math.rpi.edu/~zukerm/cgi-bin/efiles.cgi?T=37#LOOP*/
int LoopSize_energy(int** DATA, int start, int end, double* energy)
{
  double loop_energy;
  int i;
  int size_count =0;
  int plusLOOP=0; // 0:for bulge 1:for loop
  int minusLOOP=0; // 0:for bulge 1:for loop

  for (i=start+1;i<end; i++)
    {
      if (DATA[PLUS][i] != SPACE)
	{
	  plusLOOP++;
	  size_count++;
	}
    }
  for (i=start+1;i<end; i++)
    {
      if (DATA[MINUS][i] != SPACE)
	{
	  minusLOOP++;
	  size_count++;
	}
    }

  if (plusLOOP ==  minusLOOP) //loop
    {
      if (size_count <= MAX_LOOP)
	loop_energy = loop[size_count];
      else
	loop_energy = loop[MAX_LOOP] + 1.75 * R  * (273.15+37) * log((size_count+0.0)/ MAX_LOOP);
    }
  else //bulge
    {
      if (size_count <= MAX_LOOP)
	loop_energy = bulge[size_count];
      else
	loop_energy = bulge[MAX_LOOP] + 1.75 * R  * (273.15+37) * log((size_count+0.0)/ MAX_LOOP);
    }
  //debugging
  //printf("loopsize %f\n", loop_energy);
  (*energy) += loop_energy;
  return (plusLOOP && minusLOOP);  // 0:for bulge 1:for loop
}

int Min(int c, int p, int m)
{
  if ((c<=p) && (c<=m))
    return c;
  if ((p<=c) && (p<=m))
    return p;
  else
    return m;
}

/*
  asymmetry_energy = |N1 - N2| *f(min(5, N1, N2))
  N1: plus strand open loop size
  N2: minus strand open loop size
  Peritz et. al. Biochemistry 1991 30, 6428-6436
*/
double Asymmetry_energy(int** DATA,int start, int end)
{
  double energy =0;
  int i;
  int plus_size =0; 
  int minus_size=0 ; 
  int size_diff =0 ;

  for (i=start+1; i<=end-1; i++)
    if (DATA[PLUS][i] != 0)
      plus_size++;
  for (i=start+1; i<=end-1; i++)
    if (DATA[MINUS][i] != 0)
      minus_size++;
  size_diff = abs(plus_size- minus_size);

  if (size_diff)
    energy = size_diff * F_function[Min(5, plus_size, minus_size)];
  
  //debugging
  //printf("asymmetry %f\n", energy);
  return energy;
}


double Loop_Stacking (int** DATA, int start, int end, int LOOP)
{
  double energy = 0;

  if (LOOP)
    {
      /*initial mismatch stacking energy */
      energy +=  WC_Matrix[DATA[PLUS][start]][DATA[PLUS][start+1]][DATA[MINUS][start]][DATA[MINUS][start+1]]; 
      /*terminal mismatch stacking energy */
      energy +=  WC_Matrix[DATA[PLUS][end-1]][DATA[PLUS][end]][DATA[MINUS][end-1]][DATA[MINUS][end]];
    }
  else if (end-start == 2)
    /*single nucleotide bulge terminal base pair stacking energy */
    energy +=  WC_Matrix[DATA[PLUS][start]][DATA[PLUS][end]][DATA[MINUS][start]][DATA[MINUS][end]]; 
  return energy;
}


/*http://bioinfo.math.rpi.edu/~zukerm/rna/energy/node2.html*/
double Loop_Energy (int start, int end, int** DATA)
{
  double energy=0;
  int LOOP;

  if (start>= end)
    {
      printf("CODE MISTAKE\n");
      return 0;
    }
  //  printf ("%d %d \n", start, end);
  

  LOOP = LoopSize_energy(DATA,start, end, &energy);
  if ( LOOP )
    energy += Asymmetry_energy(DATA,start, end);
  energy +=Loop_Stacking (DATA, start, end, LOOP)  ;

  //debugging
  //printf("loop %f\n", energy);
  return energy;
}








