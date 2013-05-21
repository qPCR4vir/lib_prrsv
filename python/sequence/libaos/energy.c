#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "data.h"

#define MATCH 0 //middle line in blastoutput
#define PLUS 1 //first line in blastoutput
#define MINUS 2  //third line in blastoutput

#define bool int

double Perfect_Energy (int start, int end, int** DATA);
double Loop_Energy (int start, int end, int** DATA);
double Initial_Parameter();

double WC_Matrix[5][5][5][5];

int RevertCode(int code)
{
  switch(code) {
  case A:
    return T;
  case T:
    return A;
  case G:
    return C;
  case C:
    return G;
  case SPACE:
    return SPACE;
  default:
    return -1; //error
  }
}

//identify if the two code are matching base pairs eg: A:1 with T:2 and C:3 with G:4
int Matchcode (char i, char j)
{
  
  int total = i+j;
  if ( (total == (A + T) ) || (total== (G+ C)) )
    {
      if (i!= SPACE &&  j!= SPACE)
	return 1; //match
      else
	return 0; //not match
    }
  else
    return 0;  //not match
}

//code ch to int format
int CharToInt (char ch)
{
  switch (toupper(ch)) {
  case 'A':
    return A;
  case 'T':
    return T;
  case 'G':
    return G;
  case 'C':
    return C;
  case '-':
    return SPACE;
  default:
    if (isalpha(ch))
	return 'X';
    else
      return ' ';
  }
}	

/*filter input sequence to make sure every basepair is in ATGC- format
 *load two lines of sequences to data structure int* DATA[3]
 */
int LoadSeq(char* seqPlus, char* seqMinus, int** DATA, int standard_size)
{
  int i;
  int code;
  int plusSize =0;
  int minusSize =0;
  //  int* codeSeq1;
  //  int* codeSeq2;
  //  MatrixADT Matrix = (MatrixCDT *) malloc(sizeof(MatrixCDT));
  //  PathADT  max_cell = (PathCDT *)malloc(sizeof(PathCDT));

  if( (strlen(seqPlus) != standard_size) || (strlen(seqMinus) != standard_size))
    return -1; //error

  else
    {
      /*free memory allocation*/
      for (i=0; i<3; i++)
	{
	  DATA[i] = (int*) malloc(sizeof(int)*standard_size);
	  memset(DATA[i], 0 ,(standard_size));
	}

      for (i =0; i<standard_size; i++)
	{
	  //plus string
	  code = CharToInt (seqPlus[i]);
	  if (code == ' ')
	    return -1; //error
	  if (code == 'X')
	    {
	      code = CharToInt(seqMinus[i]);
	      if (code == ' ')
		return -1; //error
	      if ( (code == 'X') || (code ==SPACE))
		code = G;
	      code = RevertCode(code);
	    }
	  DATA[PLUS][i]= code;
	  if (code != SPACE)
	    ++plusSize;

	  //minus string
	  code = CharToInt (seqMinus[i]);
	  if (code == ' ')
	    return -1; //error
	  if (code == 'X')
	    {
	      code = DATA[PLUS][i];
	      if ( code == SPACE)
		code = G;
	      else
		code = RevertCode(code);
	    }
	  DATA[MINUS][i]= code;
	  if (code != SPACE)
	    ++minusSize;
	}

      return standard_size;
    }
}


void Match_Sequence(int** DATA, int standard_size)
{
  int i;

  free (DATA[MATCH]);
  DATA[MATCH] = (int*) malloc(sizeof(int)*standard_size);
  for (i =0; i<standard_size; i++)
    {
      if ( Matchcode(DATA[PLUS][i],DATA[MINUS][i]))
	DATA[MATCH][i] = '|';
      else
	DATA[MATCH][i] = ' ';
    }
}

void Print_version(int* digit_seq, int size, char* char_seq)
{
  int i;

  for (i=0; i<size; i++)
    {
      switch (digit_seq[i]) {
      case A:
	char_seq[i] = 'A';
	break;
      case T:
	char_seq[i] = 'T';
	break;
      case G:
	char_seq[i] = 'G';
	break;
      case C:
	char_seq[i] = 'C';
	break;
      case SPACE:
	char_seq[i] = '-';
	break;
      case '|':
	char_seq[i] = '|';
	break;
      case ' ':
	char_seq[i] = ' ';
	break;
      default:
	break;
      }
    }
  char_seq[size] = '\0';
}

int Revert(char* seq)
{
  int i;

  for (i=0; i<strlen(seq);i++)
    {
      switch(toupper(seq[i])) {
      case 'A':
	seq[i] = 'T';
	break;
      case 'T':
	seq[i] = 'A';
	break;
      case 'C':
	seq[i] = 'G';
	break;
      case 'G':
	seq[i] = 'C';
	break;
      case '-':
	break;
      default:
	if (isalpha(seq[i]))
	  {
	    seq[i] = 'X';
	    break;
	  }
	else
	  return 1;
      }
    }
  return 0;
}

double energy(char *s1, char *s2, int revert)
{
  /* the size of the matching pair */
  int size = strlen(s1) ;
  /* the storage int arrays for every basepair addition*/
  double addStore[size];
  double subStore[size];
  /* define the stage is perfect match = 1 or not perfect match =0 */
  //bool waston_crick=1 ;  
  
  //      int start = 0;
  //      int end =0 ;
  int loop_start =0;
  int loop_end =0;
  int loop =0;
  int i;
  double total_energy =0;
  double energy =0.0;
  int *DATA[3];


  /*if necessary, revert minusline eg ATGC-> TACG*/
  if (revert)
    if (Revert(s2))
      {
	printf ("\nWARNING: MINUSLINE SEQUENCE IS NOT IN THE CORRECT FORMAT\n");
	return 1; //error
      }

  LoadSeq(s1, s2, DATA, size);

  //match sequences
  Match_Sequence(DATA,size);

  //energy caculation
  //initialization of energy data
  total_energy = Initial_Parameter();
  //printf ("intitial %f\n", total_energy);
  
  //going down the sequence to calculate energy update
  addStore[0] = 0;
  for (i=1; i<size; i++)
    {
      if ( (DATA[MATCH][i] == '|')  && (DATA[MATCH][i-1] == '|'))
	{
	  energy = WC_Matrix[DATA[PLUS][i-1]][DATA[PLUS][i]][DATA[MINUS][i-1]][DATA[MINUS][i]];
	  addStore[i] = energy;
	  subStore[i-1] = -energy;
	}
      else if ( (DATA[MATCH][i] == ' ')  && (DATA[MATCH][i-1] == ' '))
	{
	  addStore[i] =0;
	  subStore[i-1] =0;
	}
      else if ( (DATA[MATCH][i] == '|')  && (DATA[MATCH][i-1] == ' '))
	{
	  if (i!=size-1)
	    {
	      if (DATA[MATCH][i+1] == '|')
		{
		  loop_end = i;
		  energy =Loop_Energy (loop_start, loop_end, DATA);
		  addStore[i] = energy;
		  subStore[i-1]=0;
		  subStore[loop_start] = -energy;
		  loop = 0;
		}
	      else
		{
		  addStore[i] =0;
		  subStore[i-1] =0;
		}
	    }
	  else
	    {
	      addStore[i] =0;
	      subStore[i-1] =0;
	      subStore[loop_start] = 0;
	    }
	}
      else if ( (DATA[MATCH][i] == ' ')  && (DATA[MATCH][i-1] == '|'))
	{
	  addStore[i] =0;
	  if (loop == 0)
	    {
	      loop_start = i-1; /// A-----B start : means the A :  
	      loop =1;
	    }
	  else
	    subStore[i-1] =0;
	}
      
    }
  subStore[size-1] = 0;
  
  
  for (i=0; i<size; i++)
    {
      total_energy+= addStore[i];
      //printf("%f %f\n", addStore[i], subStore[i]);
    }
  //printf("%f", total_energy);
  
  //free memeory of DATA
  for (i=0; i<3; i++)
    free(DATA[i]);
  return total_energy;
}


/* int main(int argc, char *argv[]) */
/* { */
/*   if (argc != 4) */
/*     { */
/*       printf ("\nusage: program_name first_line second_line revert(1:revert_secondline 0:no reversion)\n"); */
/*       return 1; //error */
/*     } */
/*   if (strlen(argv[1]) <2) */
/*     { */
/*       printf ("\nWARNING: SEQUENCE IS TOO SHORT\n"); */
/*       return 1; //error */
/*     } */
/*   if (*argv[3] != '1' && *argv[3] != '0' ) */
/*     { */
/*       printf ("\nWARNING: REVERT PARAMETER DID NOT SET CORRECTLY\n"); */
/*       return 1; */
/*     } */
/*   else */
/*     { */
/*       /\* */
/* 	main data structure for input data  */
/* 	data[0] is the string for matching line in blastouput */
/* 	data[1] is the string for first ouputline in blastouput */
/* 	data[2] is the string for second ouputline in blastouput */
/*       *\/ */
/*       int revert = atoi(argv[3]); */
/*       int* DATA[3]; */
/*       /\* the size of the matching pair *\/ */
/*       int size = strlen(argv[1]) ; */
/*       /\* the storage int arrays for every basepair addition*\/ */
/*       double addStore[size]; */
/*       double subStore[size]; */
/*       /\* define the stage is perfect match = 1 or not perfect match =0 *\/ */
/*       //bool waston_crick=1 ;   */
      
/*       //      int start = 0; */
/*       //      int end =0 ; */
/*       int loop_start =0; */
/*       int loop_end =0; */
/*       int loop =0; */
/*       int i; */
/*       double total_energy =0; */
/*       double energy =0.0; */
      
/*       //      char* sequence = NULL; //temporary string for holding sequence for print */

/*       /\*if necessary, revert minusline eg ATGC-> TACG*\/ */
/*       if (revert) */
/* 	if (Revert(argv[2])) */
/* 	  { */
/* 	    printf ("\nWARNING: MINUSLINE SEQUENCE IS NOT IN THE CORRECT FORMAT\n"); */
/* 	    return 1; //error */
/* 	  } */

/*       //load data into DATA (char* [3]) and filting sequences  */
/*       size = LoadSeq(argv[1], argv[2], DATA, size); */
/*       if(size == -1) */
/* 	{ */
/* 	  printf ("\nWARNING: SEQUENCE IS NOT IN THE CORRECT FORMAT\n"); */
/* 	  return 1; //error */
/* 	} */
      
/*       //match sequences */
/*       Match_Sequence(DATA,size); */

/*       //energy caculation */
/*       //initialization of energy data */
/*       total_energy = Initial_Parameter(); */
/*       //printf ("intitial %f\n", total_energy); */

/*       //going down the sequence to calculate energy update */
/*       addStore[0] = 0; */
/*       for (i=1; i<size; i++) */
/* 	{ */
/* 	  if ( (DATA[MATCH][i] == '|')  && (DATA[MATCH][i-1] == '|')) */
/* 	    { */
/* 	      energy = WC_Matrix[DATA[PLUS][i-1]][DATA[PLUS][i]][DATA[MINUS][i-1]][DATA[MINUS][i]]; */
/* 	      addStore[i] = energy; */
/* 	      subStore[i-1] = -energy; */
/* 	    } */
/* 	  else if ( (DATA[MATCH][i] == ' ')  && (DATA[MATCH][i-1] == ' ')) */
/* 	    { */
/* 	      addStore[i] =0; */
/* 	      subStore[i-1] =0; */
/* 	    } */
/* 	  else if ( (DATA[MATCH][i] == '|')  && (DATA[MATCH][i-1] == ' ')) */
/* 	    { */
/* 	      if (i!=size-1) */
/* 		{ */
/* 		  if (DATA[MATCH][i+1] == '|') */
/* 		    { */
/* 		      loop_end = i; */
/* 		      energy =Loop_Energy (loop_start, loop_end, DATA); */
/* 		      addStore[i] = energy; */
/* 		      subStore[i-1]=0; */
/* 		      subStore[loop_start] = -energy; */
/* 		      loop = 0; */
/* 		    } */
/* 		  else */
/* 		    { */
/* 		      addStore[i] =0; */
/* 		      subStore[i-1] =0; */
/* 		    } */
/* 		} */
/* 	      else */
/* 		{ */
/* 		  addStore[i] =0; */
/* 		  subStore[i-1] =0; */
/* 		  subStore[loop_start] = 0; */
/* 		} */
/* 	    } */
/* 	  else if ( (DATA[MATCH][i] == ' ')  && (DATA[MATCH][i-1] == '|')) */
/* 	    { */
/* 	      addStore[i] =0; */
/* 	      if (loop == 0) */
/* 		{ */
/* 		  loop_start = i-1; /// A-----B start : means the A :   */
/* 		  loop =1; */
/* 		} */
/* 	      else */
/* 		subStore[i-1] =0; */
/* 	    } */
	  
/* 	} */
/*       subStore[size-1] = 0; */


/*       for (i=0; i<size; i++) */
/* 	{ */
/* 	  total_energy+= addStore[i]; */
/* 	  printf("%f %f\n", addStore[i], subStore[i]); */
/* 	} */

/*       /\*old energy code starts */
/*       for (i=1; i<size;i++) */
/* 	{ */
/* 	  if (((waston_crick) && (DATA[MATCH][i] != '|')) || ((!waston_crick) && (DATA[MATCH][i] != ' ')))  */
/* 	    { */
/* 	      if (waston_crick) */
/* 		energy = Perfect_Energy (start, end, DATA); */
/* 	      else  */
/* 		{ */
/* 		  ++end; */
/* 		  if ( (size -i) < 3) */
/* 		    continue; */
/* 		  if ((DATA[MATCH][i+1]) != '|' || ((DATA[MATCH][i+2]) != '|')) */
/* 		    continue; */
/* 		  energy = Loop_Energy (start, end, DATA); */
/* 		} */
/* 	      total_energy += energy; */
/* 	      //debugging */
/* 	      //printf("update %f\n", total_energy); */
/* 	      waston_crick = !(waston_crick); */
/* 	      start = end; */
/* 	      end = i; */
/* 	    } */
/* 	  else */
/* 	    end++; */
/* 	} */

/*       //the final segment energy */
/*       if (!waston_crick) */
/* 	{ */
/* 	  //  printf("open end loop %f\n", 0.0); */
/* 	  // printf("update %f\n", total_energy); */
/* 	} */
/*       else if (end> start) */
/* 	{ */
/* 	  energy = Perfect_Energy (start, end, DATA); */
/* 	  total_energy += energy; */
/* 	  //debugging */
/* 	  //	  printf("update %f\n", total_energy); */
/* 	} */
/*       old energy code ends*\/ */

/*       /\*print */
/*       sequence = (char*) malloc(sizeof(char)*(size+1)); */
/*       memset(sequence,'\0', size+1); */
/*       printf ("\n"); */
/*       Print_version(DATA[PLUS], size,sequence); */
/*       printf ("%s\n", sequence); */

/*       Print_version(DATA[MATCH], size,sequence); */
/*       printf ("%s\n", sequence); */


/*       Print_version(DATA[MINUS], size,sequence); */
/*       printf ("%s\n", sequence); */
/*       free(sequence); */
/*       *\/ */

/*       printf("%f", total_energy); */

/*       //free memeory of DATA */
/*       for (i=0; i<3; i++) */
/* 	  free(DATA[i]); */
/*       return 1; */
/*     } */
/* } */









