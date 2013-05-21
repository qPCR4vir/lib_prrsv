/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006
Whitehead Institute for Biomedical Research, Steve Rozen
(http://jura.wi.mit.edu/rozen), and Helen Skaletsky
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

   * Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
   * Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the
distribution.
   * Neither the names of the copyright holders nor contributors may
be used to endorse or promote products derived from this software
without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "dpal.h"

/* 
 * Panic messages for when the program runs out of memory.
 */
#define DPAL_OOM_MESSAGE "Out of memory in function defined in dpal.c\n"
#define DPAL_OOM_LEN  44
#define DPAL_OOM_ERROR write(2, DPAL_OOM_MESSAGE, DPAL_OOM_LEN), exit(-2)

/*
 * We should probably remove the DPAL_FORGET_PATH compile-time option.
 * Efficiency now derives primarily for specialized versions of _dpal* for
 * particular parameter values.
 */
#ifndef DPAL_FORGET_PATH
/* 
 * Print an alignment on stderr, given the 2 aligned sequences, the "trace"
 * matrix, and the coordinate of the end point of the alignment to print.
 */
static void print_align( char *,  char *,
			int[DPAL_MAX_ALIGN][DPAL_MAX_ALIGN][3], int, int,
			 dpal_args*);
#endif

/* 
 * Return 1 if there is an illegal character in the first argument, and place
 * the illegal character in the address contained in the last argument.
 */
static int  illegal_char( char *,  dpal_ssm, char *);

static  char *xlate_ambiguity_code(char);

/*
 * The next function headers are for various versions of the
 * dynamic-programming alignment code optimized for particular input argument
 * values. 
 */

void
dpal_set_default_nt_args(dpal_args *a)
{
    unsigned int i, j;

    memset(a, 0, sizeof(*a));
    for (i = 0; i <= UCHAR_MAX; i++)
	for (j = 0; j <= UCHAR_MAX; j++)
	    if (('A' == i || 'C' == i || 'G' == i || 'T' == i || 'N' == i)
		&& ('A' == j || 'C' == j || 'G' == j || 'T' == j 
		    || 'N' == j)) {
		    if (i == 'N' || j == 'N') 
			a->ssm[i][j] = -25;
		    else if (i == j)
			a->ssm[i][j] = 100;
		    else 
			a->ssm[i][j] = -100;
		} else
		    a->ssm[i][j] = INT_MIN;

    a->check_chars        = 1;
    a->debug              = 0;
    a->fail_stop          = 1;
    a->flag               = DPAL_LOCAL;
    a->force_generic      = 0;
    a->force_long_generic = 0;
    a->force_long_maxgap1 = 0;
    a->gap                = -100;
    a->gapl               = -100;
    a->max_gap            = 3;
    a->score_only         = 0;
}

void
dpal_set_h_nt_matrix(dpal_args *a)
{
    unsigned int i, j;

    for (i = 0; i <= UCHAR_MAX; i++)
	for (j = 0; j <= UCHAR_MAX; j++)
	    if (('A' == i || 'C' == i || 'G' == i || 'T' == i || 'N' == i)
		&& ('A' == j || 'C' == j || 'G' == j || 'T' == j 
		    || 'N' == j)) {
		    if (i == 'N' || j == 'N') 
			a->ssm[i][j] = -50;
		    else if (i == j) {
		      if ('C' == i || 'G' == i)
			a->ssm[i][j] = 300;
		      else
			a->ssm[i][j] = 200;
		    }
		    else 
			a->ssm[i][j] = -50;
		} else
		    a->ssm[i][j] = INT_MIN;
}

/* The argument a must be a DNA scoring matrix.
   Modify a so that it for a match between
   any two ambiguity codes (or between ambiguity code and base),
   e.g. B and S, the score will be the maximum of
   score between any base in B and any base in S,
   in the example between any pair in {C, G, T} X {C, G}.
*/
int
dpal_set_ambiguity_code_matrix(dpal_args *a)
{
   char *c1, *c2;
   char *amb_codes = "BDHVRYKMSWN";
   char *all_bases ="ACGT";
   char *bases1, *bases2, *b1, *b2;
  int extreme;
  for (c1 = amb_codes; *c1; c1++) {
    bases1 = xlate_ambiguity_code(*c1);
    if (!bases1) return 0;

    /* Do matches between c1 and all other
       ambiguity codes. */
    for (c2 = amb_codes; *c2; c2++) {
      bases2 = xlate_ambiguity_code(*c2);
      if (!bases2) return 0;
      extreme = INT_MIN;
      for (b1 = bases1; *b1; b1++) {
	for (b2 = bases2; *b2; b2++) {
	  if (a->ssm[*b1][*b2] > extreme) {
	    extreme = a->ssm[*b1][*b2];
	  }
	}
      }
      /* extreme is now the maximum score
	 for a match between any 2 bases
         represented *c1, *c2. */
      a->ssm[*c1][*c2] = extreme;
    }

    /* Do matches between c1 and all bases. */
    for (b2 = all_bases; *b2; b2++) {
      extreme = INT_MIN;
      for (b1 = bases1; *b1; b1++) {
	if (a->ssm[*b1][*b2] > extreme) {
	  extreme = a->ssm[*b1][*b2];
	}
      }
      a->ssm[*c1][*b2] = extreme;
      a->ssm[*b2][*c1] = extreme;
    }
  }
  return 1;
}

static  char *
xlate_ambiguity_code(char c)
{
  if ('N' == c)      return "ACGT";
  else if ('B' == c) return "CGT";
  else if ('D' == c) return "AGT";
  else if ('H' == c) return "ACT";
  else if ('V' == c) return "ACG";
  else if ('R' == c) return "AG";
  else if ('Y' == c) return "CT";
  else if ('K' == c) return "GT";
  else if ('M' == c) return "AC";
  else if ('S' == c) return "CG";
  else if ('W' == c) return "AT";
  else return NULL; /* Error condition */
}

#define CHECK_ERROR(COND,MSG) if (COND) { out->msg = MSG; goto FAIL; }

void
_dpal_long_nopath_maxgap1_local(
     char *X, char *Y,
     int xlen, int ylen,
     dpal_args *in,
    dpal_results *out)
{
    /* The "score matrix" (matrix of best scores). */
    int *S0, *S1, *S2; 
    int *P0, *P1, *P2;
    int *S;

    register int i, j;
    register int gap = in->gap;
    register int smax;           /* The optimum score. */
    register int score;          /* Current score. */
    register int a;

#ifdef DPAL_PRINT_COVERAGE
    fprintf(stderr, "_dpal_long_nopath_maxgap1_local called\n");
#endif

    P0 = (int*)malloc(sizeof(int)*ylen);
    if (!P0) { DPAL_OOM_ERROR; }
    P1 = (int*)malloc(sizeof(int)*ylen);
    if (!P1) { DPAL_OOM_ERROR; }
    P2 = (int*)malloc(sizeof(int)*ylen);
    if (!P2) { DPAL_OOM_ERROR; }

    S0 = P0; S1 = P1; S2 = P2;

    smax = 0; /* For local alignment score can never be less than 0. */

    /* Initialize the 0th row of the score matrix. */
    for(j=0; j < ylen; j++) { 
	score = in->ssm[X[0]][Y[j]]; 
	if (score < 0) score = 0;
	else if (score > smax) smax = score;
	/*S[0][j] = score;*/
	S0[j] = score;
    }	

    /* Set the 1st row of the score matrix. */
    score = in->ssm[X[1]][Y[0]];
    if(score < 0) score = 0;
    else if (score > smax) smax = score;
    S1[0] = score;
    for(j=1; j < ylen; j++) {
	score = S0[j-1];
	if(j>1 && (a=S0[j-2] + gap) > score)score = a;
	score += in->ssm[X[1]][Y[j]];
	if (score < 0) score = 0;
	else if(score > smax) smax = score;
	S1[j] = score;
    }

    for(i=2; i < xlen; i++) {
	score = in->ssm[X[i]][Y[0]];
	if (score < 0) score = 0;
	else if (score > smax) smax = score;
	S2[0] = score;
	score = S1[0];
	if((a=S0[0] + gap) > score) score = a;
	score += in->ssm[X[i]][Y[1]];
	if(score < 0) score = 0;
	else if (score > smax) smax = score;
	S2[1] = score;
	for(j=2; j < ylen; j++) {
	    score = S0[j-1];
	    if((a=S1[j-2])>score) score = a;
	    score +=gap;
	    if((a=S1[j-1]) >score) score = a;

	    score += in->ssm[X[i]][Y[j]];	
	    if (score < 0 ) score = 0;
	    else if (score > smax) smax = score;
	    S2[j]=score;
	}
	S = S0; S0 = S1; S1 = S2; S2 = S;
    }
    out->score = smax;
    out->path_length=0;
    free(P0); free(P1); free(P2);
} /* _dpal_long_nopath_maxgap1_local */

void
_dpal_long_nopath_maxgap1_global_end(
     char *X, char *Y,
     int xlen, int ylen,
     dpal_args *in,
    dpal_results *out)
{
    /* The "score matrix" (matrix of best scores). */
    int *S0, *S1, *S2, *S; 
    int *P0, *P1, *P2;

    register int i, j, k;
    register int gap = in->gap;
    register int smax;           /* The optimum score. */
    register int score;          /* Current score. */
    register int a, t;

#ifdef DPAL_PRINT_COVERAGE
    fprintf(stderr, "_dpal_long_nopath_maxgap1_global_end called\n");
#endif

    P0 = (int*)malloc(sizeof(int)*xlen);
    if (!P0) { DPAL_OOM_ERROR; }
    P1 = (int*)malloc(sizeof(int)*xlen);
    if (!P1) { DPAL_OOM_ERROR; }
    P2 = (int*)malloc(sizeof(int)*xlen);
    if (!P2) { DPAL_OOM_ERROR; }

    S0 = P0; S1 = P1; S2 = P2;

    smax = in->ssm[X[xlen-1]][Y[0]];
	   
    /* Set the 0th row of the score matrix. */
    for(j=0; j<xlen; j++) S0[j] = in->ssm[X[j]][Y[0]]; 

    /* Set the 1st row of the score matrix. */
    S1[0] = in->ssm[X[0]][Y[1]];
    for(j=1; j < xlen; j++){
	score = S0[j-1];
	if(j>1 && (a=S0[j-2] + gap)> score)score = a;
	score += in->ssm[X[j]][Y[1]];
	if(score > smax && j == xlen-1) smax = score;
        S1[j] = score;
    }

    k = ylen - (int)(xlen / 2) + 1;
    if (k<1) k = 1;

    /* Set the rectangular part of almost the remainder of the matrix. */
    for(j=2; j<k+1; j++) {
	S2[0] = in->ssm[X[0]][Y[j]];
	score = S1[0];
	if((a=S0[0]+gap) > score) score = a;
	score += in->ssm[X[1]][Y[j]];
	S2[1] = score;
       for(i=2; i<xlen-1; i++) {
	   score = S1[i-2];
	   if((a=S0[i-1]) > score)score = a;
	   score += gap;
	   if((a=S1[i-1]) > score)score = a;
	   score += in->ssm[X[i]][Y[j]];
	   S2[i] = score;
        }
	score = S1[xlen-3];
	if((a=S0[xlen-2]) > score)score = a;
	score += gap;
	if((a=S1[xlen-2]) > score)score = a;
	score += in->ssm[X[xlen-1]][Y[j]];
	S2[xlen-1] = score;
	if(score > smax) smax = score;
	S = S0; S0 = S1; S1 = S2; S2 = S;
    }

    /* Set the triangular part of almost the remainder of the matrix. */
    t = 2;
    for(j=k+1; j<ylen; j++) {
       for(i=t; i<xlen-1; i++) {
	    score = S1[i-2];
	    if((a=S0[i-1]) > score) score = a;
	    score += gap;
	    if((a=S1[i-1]) > score) score = a;
	    score += in->ssm[X[i]][Y[j]];
	    S2[i] = score;
       }
       t += 2;
       score = S1[xlen-3];
       if((a=S0[xlen-2]) > score)score = a;
       score += gap;
       if((a=S1[xlen-2]) > score)score = a;
       score += in->ssm[X[xlen-1]][Y[j]];
       S2[xlen-1] = score;
       if(score > smax) smax = score;
       S = S0; S0 = S1; S1 = S2; S2 = S;
    }

    free(P0); free(P1); free(P2);
    out->score = smax;
    out->path_length=0;
} /* _dpal_long_nopath_maxgap_global_end */



static int
illegal_char(
     char *X,
     dpal_ssm ssm,
    char *out)
{
    register  char *p;
    for (p = X; *p != '\0' && ssm[*p][*p] != INT_MIN; p++);
    if (*p == '\0')
	return 0;
    else {
	*out = *p;
	return 1;
    }
}

