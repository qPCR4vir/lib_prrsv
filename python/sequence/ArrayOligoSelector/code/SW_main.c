#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <limits.h>

#include "SmithWaterman.h"

int Smith_Waterman_Fill (int* seq1Arr, int* seq2Arr, MatrixADT Matrix, PathADT max_cell, Score_Matrix* score_Matrix);
Score_Matrix *choose_matrix(int matrix_name);
void Output_Best_SW_Alignment (char *seq1, char *seq2, MatrixADT Matrix, PathADT max_cell, char *buff);
void free_MatrixADT(MatrixADT Matrix);
void free_PathADT (PathADT path);
void  free_score_Matrix(Score_Matrix* score_Matrix);
int* Encode_Seq (char* seq);
void Complement_Seq_RDNA (char *seq, char *comp_seq);
void InputFileToString(FILE* fin, char** seq);

int SW (char *seq1, char *seq2, int alignment, int size,char* buff)
{
  MatrixADT Matrix = (MatrixCDT *) malloc(sizeof(MatrixCDT));
  PathADT  max_cell;
  Score_Matrix *score_Matrix;
  int* codeSeq1;
  int* codeSeq2;
  int max_score;
  int temp_seq[size];
  char seq[size+1];
  int i,j;
  int smatrix_name = PAM47;

  /* score matrix */
  score_Matrix = choose_matrix(smatrix_name);
  if (score_Matrix == NULL)
    {
      free_score_Matrix(score_Matrix);
      score_Matrix= NULL;
      return;
    }

  Matrix->row = size +1;
  Matrix->column = size+1;
  Matrix->grid = (Matrix_Cell ***) malloc(sizeof(Matrix_Cell **)*Matrix->row);
  for (i=0; i< Matrix->row; i++)
    {
      Matrix->grid[i] = (Matrix_Cell **) malloc(sizeof(Matrix_Cell *)*Matrix->column);
    }

  for (i=0; i< Matrix->row; i++)
    {
      for (j=0; j< Matrix->column; j++)
	{
	  Matrix->grid[i][j] = (Matrix_Cell *) malloc (sizeof(Matrix_Cell));
	}
    }

  //encode seq1 and seq2
  codeSeq1 = Encode_Seq(seq1);
  codeSeq2 = Encode_Seq(seq2);

  for (i =0; i<strlen(seq1)-size+1; i++)
    {
      //reverse seq2 segment
      for (j =0; j<size; j++)
	{
	  temp_seq[j] = (codeSeq2+i) [size-j-1];
	}
      max_cell = (PathCDT *)malloc(sizeof(PathCDT));
      max_cell->start = NULL;
      max_score = Smith_Waterman_Fill (codeSeq1+i, temp_seq, Matrix, max_cell,score_Matrix);

      if ((alignment) && (strlen(seq2) == size))
	{
	  for (i=0; i<size; i++)
	    seq[i] = seq2[size-i-1];
	  seq[size] = '\0';
	  Output_Best_SW_Alignment(seq1, seq, Matrix, max_cell,buff);
	}

      free_PathADT(max_cell);
      sprintf (buff+strlen(buff),"%d\n", max_score);
    }
  

  free(codeSeq1);
  free(codeSeq2);
  
  free_MatrixADT(Matrix);

  //free score_Matrix memory
  free_score_Matrix(score_Matrix);
  score_Matrix = NULL;

  return max_score;
}

int main (int argc, char *argv[])
{
  FILE *fin;
  char *seq ;
  char *comp_seq;
  int ALIGNMENT;
  int size;
  char *buff;
  int i;

  if (argc !=  4)
    {
      printf ("\nusage: excutable seqstringfile print_alignment(0:1) segment_size\n");
      return 0; //error
    }
  
  else
    {
      //seq = argv[1];
      fin = fopen(argv[1],"r");
      if (fin == NULL)
	return 0;
      InputFileToString(fin, &seq) ;
          
      ALIGNMENT = atoi(argv[2]);
      size = atoi(argv[3]);
      
      comp_seq = (char *)malloc(sizeof(char)*(strlen(seq)+1));
      memset(comp_seq, 0, (strlen(seq)+1));
      
      buff=(char *)malloc(sizeof(char)*(strlen(seq)+1)*5);
      memset(buff, 0, (strlen(seq)+1));
      //printf("%d\n",strlen(buff));

      //complement sequence not reversed yet!
      Complement_Seq_RDNA (seq, comp_seq); 
      if (comp_seq == NULL)
	{
	  free (comp_seq);
	  return 0; //error
	}

      SW (seq, comp_seq, ALIGNMENT,size,buff);
      //printf("%d\n",strlen(buff));
      
      printf("%s",buff);
      //write(stdout, buff , strlen( buff ));
      //printf("%d\n",strlen(buff));
      //free comp_seq memory
      free (comp_seq);
      free(seq);

      return 1; //success
    }
}
