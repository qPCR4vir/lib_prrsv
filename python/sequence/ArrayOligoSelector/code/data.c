#include <stdio.h>
#include "data.h"

double WC_Matrix[5][5][5][5];
double initial_WC_pair = 3.4;
double loop[MAX_LOOP+1];
double bulge[MAX_LOOP+1];

int i, j, k, m;

/* All energy parameters are in kcal/mol unit 
   All energy parameters are in 37 degree C standard condition
*/ 


/*N. Sugimoto et.al.  NAR  1996, Vol 24, No22 4501-4505*/
void Initial_WC_Parameter()
{
  for(i=0; i<5; i++)
    {
      for(j=0; j<5; j++)
	{
	  for( k=0; k<5; k++)
	    {
	      for (m=0; m<5; m++)
		{
		  WC_Matrix[i][j][k][m] = 0.0;
		}
	    }
	}
    }
  
  WC_Matrix[A][A][T][T] = -1.2;
  WC_Matrix[T][T][A][A] = -1.2;

  WC_Matrix[A][T][T][A] = -0.9;

  WC_Matrix[T][A][A][T] = -0.9;

  WC_Matrix[C][A][G][T] = -1.7;
  WC_Matrix[T][G][A][C] = -1.7;

  WC_Matrix[C][T][G][A] = -1.5;
  WC_Matrix[A][G][T][C] = -1.5;

  WC_Matrix[G][A][C][T] = -1.5;
  WC_Matrix[T][C][A][G] = -1.5;

  WC_Matrix[G][T][C][A] = -1.5;
  WC_Matrix[A][C][T][G] = -1.5;

  WC_Matrix[C][G][G][C] = -2.8;

  WC_Matrix[G][C][C][G] = -2.3;

  WC_Matrix[G][G][C][C] = -2.1;
  WC_Matrix[C][C][G][G] = -2.1;
}

/*N. Peyret et.al. Biochemistry 1999, 38, 3468-3477
  Biochemistry Vol37, No 26, 1998
  NAR 1998, Vol 26, No.11
  Biochemistry, Vol 37, No. 8, 1998
  Biochemistry, Vol 36, No. 34, 1997
*/
void Initial_mismatch_Parameter()
{
  WC_Matrix[A][A][T][A] = 0.61;
  WC_Matrix[A][T][A][A] = 0.61;
  WC_Matrix[C][A][G][A] = 0.43;
  WC_Matrix[A][G][A][C] = 0.43;
  WC_Matrix[G][A][C][A] = 0.17;
  WC_Matrix[A][C][A][G] = 0.17;
  WC_Matrix[T][A][A][A] = 0.69;
  WC_Matrix[A][A][A][T] = 0.69;

  WC_Matrix[A][C][T][C] = 1.33;
  WC_Matrix[C][T][C][A] = 1.33;
  WC_Matrix[C][C][G][C] = 0.70;
  WC_Matrix[C][G][C][C] = 0.70;
  WC_Matrix[G][C][C][C] = 0.79;
  WC_Matrix[C][C][C][G] = 0.79;
  WC_Matrix[T][C][A][C] = 1.05;
  WC_Matrix[C][A][C][T] = 1.05;

  WC_Matrix[A][G][T][G] = -0.13;
  WC_Matrix[G][T][G][A] = -0.13;
  WC_Matrix[G][T][G][A] = -0.13;
  WC_Matrix[A][G][T][G] = -0.13;
  WC_Matrix[C][G][G][G] = -0.11;
  WC_Matrix[G][G][G][C] = -0.11;
  WC_Matrix[G][G][C][G] = -1.11;
  WC_Matrix[G][C][G][G] = -1.11;
  WC_Matrix[T][G][A][G] = 0.44;
  WC_Matrix[G][A][G][T] = 0.44;

  WC_Matrix[A][T][T][T] = 0.69;
  WC_Matrix[T][T][T][A] = 0.69;
  WC_Matrix[C][T][G][T] = -0.12;
  WC_Matrix[T][G][T][C] = -0.12;
  WC_Matrix[G][T][C][T] = 0.45;
  WC_Matrix[T][C][T][G] = 0.45;
  WC_Matrix[T][T][A][T] = 0.68;
  WC_Matrix[T][A][T][T] = 0.68;


  WC_Matrix[A][A][T][C] = 0.88;
  WC_Matrix[C][T][A][A] = 0.88;
  WC_Matrix[A][C][T][A] = 0.77;
  WC_Matrix[A][T][C][A] = 0.77;

  WC_Matrix[C][A][G][C] = 0.75;
  WC_Matrix[C][G][A][C] = 0.75;
  WC_Matrix[C][C][G][A] = 0.79;
  WC_Matrix[A][G][C][C] = 0.79;

  WC_Matrix[G][A][C][C] = 0.81;
  WC_Matrix[C][C][A][G] = 0.81;
  WC_Matrix[G][C][C][A] = 0.47;
  WC_Matrix[A][C][C][G] = 0.47;

  WC_Matrix[T][A][A][C] = 0.92;
  WC_Matrix[C][A][A][T] = 0.92;
  WC_Matrix[T][C][A][A] = 1.33;
  WC_Matrix[A][A][C][T] = 1.33;


  WC_Matrix[A][C][T][T] = 0.64;
  WC_Matrix[T][T][C][A] = 0.64;
  WC_Matrix[A][T][T][C] = 0.73;
  WC_Matrix[C][T][T][A] = 0.73;

  WC_Matrix[C][C][G][T] = 0.62;
  WC_Matrix[T][G][C][C] = 0.62;
  WC_Matrix[C][T][G][C] = 0.40;
  WC_Matrix[C][G][T][C] = 0.40;

  WC_Matrix[G][C][C][T] = 0.62;
  WC_Matrix[T][C][C][G] = 0.62;
  WC_Matrix[G][T][C][C] = 0.98;
  WC_Matrix[C][C][T][G] = 0.98;

  WC_Matrix[T][C][A][T] = 0.97;
  WC_Matrix[T][A][C][T] = 0.97;
  WC_Matrix[T][T][A][C] = 0.75;
  WC_Matrix[C][A][T][T] = 0.75;


  WC_Matrix[A][A][T][G] = 0.14;
  WC_Matrix[G][T][A][A] = 0.14;
  WC_Matrix[A][G][T][A] = 0.02;
  WC_Matrix[A][T][G][A] = 0.02;

  WC_Matrix[C][A][G][G] = 0.03;
  WC_Matrix[G][G][A][C] = 0.03;
  WC_Matrix[C][G][G][A] = 0.11;
  WC_Matrix[A][G][G][C] = 0.11;

  WC_Matrix[G][A][C][G] = -0.25;
  WC_Matrix[G][C][A][G] = -0.25;
  WC_Matrix[G][G][C][A] = -0.52;
  WC_Matrix[A][C][G][G] = -0.52;

  WC_Matrix[T][A][A][G] = 0.42;
  WC_Matrix[G][A][A][T] = 0.42;
  WC_Matrix[T][G][A][A] = 0.74;
  WC_Matrix[A][A][G][T] = 0.74;


  WC_Matrix[A][G][T][T] = 0.71;
  WC_Matrix[T][T][G][A] = 0.71;
  WC_Matrix[A][T][T][G] = 0.07;
  WC_Matrix[G][T][T][A] = 0.07;

  WC_Matrix[C][G][G][T] = -0.47;
  WC_Matrix[T][G][G][C] = -0.47;
  WC_Matrix[C][T][G][G] = -0.32;
  WC_Matrix[G][G][T][C] = -0.32;

  WC_Matrix[G][G][C][T] = 0.08;
  WC_Matrix[T][C][G][G] = 0.08;
  WC_Matrix[G][T][C][G] = -0.59;
  WC_Matrix[G][C][T][G] = -0.59;

  WC_Matrix[T][G][A][T] = 0.43;
  WC_Matrix[T][A][G][T] = 0.43;
  WC_Matrix[T][T][A][G] = 0.34;
  WC_Matrix[G][A][T][T] = 0.34;
}

/*http://bioinfo.math.rpi.edu/~zukerm/cgi-bin/efiles.cgi?T=37#LOOP*/
void Initial_loop_Parameter()
{
  loop[0]=0;
  loop[1]=0;
  loop[2]=             4.10;
  loop[3]=              5.10;
  loop[4]=              4.90;
  loop[5]=              5.30;
  loop[6]=              5.70;
  loop[7]=              5.90;
  loop[8]=              6.00;
  loop[9]=              6.10;
  loop[10]=              6.30;
  loop[11]=              6.40;
  loop[12]=              6.40;
  loop[13]=              6.50;
  loop[14]=              6.60;
  loop[15]=              6.70;
  loop[16]=              6.80;
  loop[17]=              6.80;
  loop[18]=              6.90;
  loop[19]=              6.90;
  loop[20]=              7.00;
  loop[21]=              7.10;
  loop[22]=              7.10;
  loop[23]=              7.10;
  loop[24]=              7.20;
  loop[25]=              7.20;
  loop[26]=              7.30;
  loop[27]=              7.30;
  loop[28]=              7.40;
  loop[29]=              7.40;
  loop[30]=              7.40;
}

/*http://bioinfo.math.rpi.edu/~zukerm/cgi-bin/efiles.cgi?T=37#LOOP*/
void Initial_bulge_Parameter()
{
  bulge[0]=0;
  bulge[1]= 3.90;
  bulge[2]=              3.10;
  bulge[3]=              3.50;
  bulge[4]=              4.20;
  bulge[5]=              4.80;
  bulge[6]=              5.00;
  bulge[7]=              5.20;
  bulge[8]=              5.30;
  bulge[9]=              5.40;
  bulge[10]=              5.50;
  bulge[11]=              5.70;
  bulge[12]=              5.70;
  bulge[13]=              5.80;
  bulge[14]=              5.90;
  bulge[15]=              6.00;
  bulge[16]=              6.10;
  bulge[17]=              6.10;
  bulge[18]=              6.20;
  bulge[19]=              6.20;
  bulge[20]=              6.30;
  bulge[21]=              6.30;
  bulge[22]=              6.40;
  bulge[23]=              6.40;
  bulge[24]=              6.50;
  bulge[25]=              6.50;
  bulge[26]=              6.50;
  bulge[27]=              6.60;
  bulge[28]=              6.70;
  bulge[29]=              6.70;
  bulge[30]=              6.70;
}


/* 
   Peritz et. al. Biochemistry 1991 30, 6428-6436
*/
double F_function[6]= {0,0.7, 0.6, 0.4, 0.1};  


double Initial_Parameter()
{
  double energy = initial_WC_pair;

  Initial_WC_Parameter();
  Initial_mismatch_Parameter();
  Initial_loop_Parameter();
  Initial_bulge_Parameter();
  return energy;
}






