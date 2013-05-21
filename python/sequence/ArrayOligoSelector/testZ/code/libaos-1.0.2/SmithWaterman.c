#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <limits.h>
#include "data.h"
#include "SmithWaterman.h"

void free_Node (Node* node);

Score_Matrix *choose_matrix(int matrix_name)
{
  Score_Matrix *MatrixPtr= NULL;
  int i, j;
  int match=0;
  int mis_match=0;

  switch (matrix_name)
    {
    case (PAM47):
      MatrixPtr = (Score_Matrix *) malloc(sizeof(Score_Matrix));
      MatrixPtr->open = -7;
      MatrixPtr->extend = 0;
      MatrixPtr->row = DNA_ALPHABET;
      MatrixPtr->column =DNA_ALPHABET;
      MatrixPtr->grid = (int **) malloc(sizeof(int *)*DNA_ALPHABET);
      match = 5;
      mis_match = -4;

      for (i = 0; i<MatrixPtr->row; i++)
	{
	  MatrixPtr->grid[i] = (int *) malloc(sizeof(int)*DNA_ALPHABET);
	}

      for (i=0; i<MatrixPtr->row; i++)
	{
	  for (j=0; j<MatrixPtr->column; j++)
	    MatrixPtr->grid[i][j] = mis_match;
	}
      MatrixPtr->grid[A][A] = MatrixPtr->grid[T][T] = MatrixPtr->grid[C][C] = MatrixPtr->grid[G][G] = match;
      MatrixPtr->grid[A][N] = MatrixPtr->grid[T][N] = MatrixPtr->grid[C][N] = MatrixPtr->grid[G][N] = match;
      MatrixPtr->grid[N][A] = MatrixPtr->grid[N][T] = MatrixPtr->grid[N][C] = MatrixPtr->grid[N][G] = match;
      MatrixPtr->grid[N][N]  = match;
      break;
    default:
      break;
    }
  return MatrixPtr;
}

void  free_score_Matrix(Score_Matrix* score_Matrix)
{
  int i;
  for (i = 0; i<score_Matrix->row; i++)
    {
      free (score_Matrix->grid[i]);
      score_Matrix->grid[i] = NULL;
    }
  free( score_Matrix->grid);
  score_Matrix->grid = NULL;
  free( score_Matrix);
  score_Matrix= NULL;
}

static int Encode (char ch)
{
  switch (toupper(ch))
    {
    case 'A':
      return A;
    case 'T':
      return T;
    case 'G':
      return G;
    case 'C':
      return C;
    case 'a':
      return A;
    case 't':
      return T;
    case 'g':
      return G;
    case 'c':
      return C;
    case 'N':
      return N;
    case 'n':
      return N;
    default:
      return N;
    }
}

/*best local alignment segment in codes*/
int Best_SW_Align_Code (PathADT max_cell, MatrixADT Matrix, int** line1, int** line2,  int* code1, int* code2)
{
  Matrix_Cell* cell;
  Node* end_node;
  int number=0, i, end_row, end_column;
  int* rowVector;
  int* columnVector;

  end_node = max_cell->start;
  if (end_node)
    {
      end_row= end_node->row;
      end_column = end_node->column;
    }
  else
      return 0;

  //size of the alignment
  cell = Matrix->grid[end_row][end_column] ;
  while ( cell->prev != NULL)
    {
      cell = cell->prev;
      number++;
    }
  
  rowVector= (int*)malloc(sizeof(int)*(number));
  columnVector = (int*)malloc(sizeof(int)*(number));
  cell = Matrix->grid[end_row][end_column] ;
  for (i=number-1; i>=0; i--)
    {
      rowVector[i] = cell->row;
      columnVector[i] = cell->column;
      cell = cell->prev;
    }

  *line1 = (int*)malloc(sizeof(int)*(number));
  *line2 = (int*)malloc(sizeof(int)*(number));
  for (i=0; i<number; i++)
    {
      if (i==0)
	{
	  (*line1)[i] = code1[rowVector[i] -1 ]; 
	  (*line2)[i] = code2[columnVector[i] -1 ];      
	}
      else
	{
	  if (rowVector[i] == rowVector[i-1])
	    (*line1)[i] = SPACE;
	  else
	    (*line1)[i] = code1[rowVector[i] -1]; 
	  if (columnVector[i] == columnVector[i-1])
	    (*line2)[i] = SPACE;
	  else
	    (*line2)[i] = code2[ columnVector[i] -1 ];      
	}
    }
  free(rowVector);
  rowVector= NULL;
  free(columnVector);
  columnVector = NULL;
  return number;
}


/*best local alignment print  */
void Output_Best_SW_Alignment (char *seq1, char *seq2, MatrixADT Matrix, PathADT max_cell, char *buff)
{
  Node *end_node;
  Matrix_Cell *start_cell, *cell;
  Matrix_Cell *revArr_cell[(strlen(seq1)+strlen(seq2))];

  int end_row=0, end_column=0, start_row =0, start_column =0, anchor=0;
  int i;
  int number=0;

  end_node = max_cell->start;

  /*alignmnet end position  if exists*/
  if (end_node)
    {
      end_row= end_node->row;
      end_column = end_node->column;
    }
  else
    return;

  /*length of the alignment*/
  //trick: the padding row or column is ignored on the counting
  start_cell = cell = Matrix->grid[end_row][end_column] ;
  while ( cell->prev != NULL)
    {
      cell = cell->prev;
      number++;
    }


  /*set the matrix_cell ptr in the alignment into an array */
  memset(revArr_cell,'\0', (strlen(seq1)+ strlen(seq2)));
  for (i=0; i<number; i++)
    {
      revArr_cell[i] = start_cell;
      start_cell = start_cell->prev;
    }

  if (number)
    {
      start_row = revArr_cell[number-1]->row -1;
      start_column = revArr_cell[number-1]->column-1;
    }

  else
    {
      sprintf(buff+strlen(buff),"\n");
      sprintf(buff+strlen(buff),"\n");
      return;
    }

  anchor = Max(start_row, start_column);

  /*space before alignment start for seq1  */
  for (i=0; i<anchor-start_row ; i++)
    sprintf(buff+strlen(buff),"%c", (' '));

  /*print seq1 sequence before anchor position */
  for (i=0; i<start_row ; i++)
    sprintf(buff+strlen(buff),"%c", ('.'));
  
  if (number)
    {
      sprintf(buff+strlen(buff),"%c",seq1[revArr_cell[number-1]->row -1]);
      for ( i=number-2 ; i>=0; i--)
	{
	  if ( revArr_cell[i]->row != revArr_cell[i+1]->row)
	    {
	      if ((seq1[revArr_cell[i]->row-1] != seq2[revArr_cell[i]->column-1]) || (revArr_cell[i]->column == revArr_cell[i+1]->column))
		sprintf(buff+strlen(buff),".");
	      else
		sprintf(buff+strlen(buff),"%c", (seq1[revArr_cell[i]->row -1]));
	    }
	  else
	    sprintf(buff+strlen(buff)," ");
	}
    }

  /*print seq1 sequence after alignment*/
  for (i=end_row; i<strlen(seq1) ; i++)
    sprintf(buff+strlen(buff),"%c", ('.'));
    
  sprintf(buff+strlen(buff),"\n");



  /*space before alignment start for seq2 */
  for (i=0; i<anchor-start_column ; i++)
    sprintf(buff+strlen(buff),"%c", (' '));

  /*print seq2 sequence before anchor position */
  for (i=0; i<start_column ; i++)
    sprintf(buff+strlen(buff),"%c", ('.'));

  if (number)
    {
      sprintf(buff+strlen(buff),"%c",seq2[revArr_cell[number-1]->column -1]);
      for ( i=number-2 ; i>=0; i--)
	{
	  if ( revArr_cell[i]->column != revArr_cell[i+1]->column)
	    {
	      if ((seq1[revArr_cell[i]->row-1] != seq2[revArr_cell[i]->column-1]) || (revArr_cell[i]->row == revArr_cell[i+1]->row))
		sprintf(buff+strlen(buff),".");
	      else
		sprintf(buff+strlen(buff),"%c", (seq2[revArr_cell[i]->column -1]));
	    }
	  else
	    sprintf(buff+strlen(buff)," ");
	}
    }
  for (i=end_column; i<strlen(seq2) ; i++)
    sprintf(buff+strlen(buff),"%c", ('.'));
    
  sprintf(buff+strlen(buff),"\n");
  return;
}

int* Encode_Seq (char* seq)
{
  int i;
  int size = strlen(seq);

  int* codeSeq = (int*) malloc(sizeof(int)* size);
  for (i=0; i<size; i++)
    {
      codeSeq[i] = Encode(seq[i]);
      if( codeSeq[i] == UNDEFINED)
	return NULL;
    }
  return codeSeq;
}

/*row in matrix seq1Arr column in matrix is seq2Arr*/
int Smith_Waterman_Fill (int* seq1Arr, int* seq2Arr, MatrixADT Matrix, PathADT max_cell, Score_Matrix *score_Matrix)
{
  int m, k ;
  int score_dia, score_ver, score_hor, cell_score, max_score=0;
  Node *temp_node, *node;
  int i, j;
  int open_gap , extention_gap;


  m= Matrix->row; 
  k = Matrix->column;   /* m : row   n: column */

  /*Matrix initialization */
  //padding
  for (i=0; i< m; i++)
    {
      Matrix->grid[i][0]->score = 0;
      Matrix->grid[i][0]->gap =NOGAP;
      Matrix->grid[i][0]->prev = NULL;
      Matrix->grid[i][0]->row =i;
      Matrix->grid[i][0]->column = 0;
    }
  //padding
  for (i=0; i< k; i++)
    {
      Matrix->grid[0][i]->score = 0;
      Matrix->grid[0][i]->gap =NOGAP;
      Matrix->grid[0][i]->prev = NULL;
      Matrix->grid[0][i]->row = 0;
      Matrix->grid[0][i]->column = i;
    }
  
  /*left to right, up to down fill Matrix */
  for (i =1; i<m; i ++) //row in matrix
    {
      for (j=1; j< k; j++) //column in matrix
	{
	  score_dia = Matrix->grid[i-1][j-1]->score + score_Matrix->grid[seq1Arr[i-1]] [seq2Arr[j-1]];

	  open_gap =0;
	  extention_gap =0;
          if ( Matrix->grid[i-1][j]->gap == V_GAP)
	    extention_gap =1;
	  else if (Matrix->grid[i-1][j]->score != 0)
	    open_gap =1;
	  score_ver = Matrix->grid[i-1][j]->score + extention_gap*(score_Matrix->extend)  + open_gap * ( score_Matrix->open); 

	  open_gap =0;
	  extention_gap =0;
	  if ( Matrix->grid[i][j-1]->gap == H_GAP)
	    extention_gap =1;
	  else if (Matrix->grid[i][j-1]->score != 0)
	    open_gap =1;
	  score_hor =Matrix->grid[i][j-1]->score + extention_gap*score_Matrix->extend + open_gap*score_Matrix->open;

	  /* Maxi */
	  cell_score = Max (score_dia, score_ver);
	  cell_score = Max (cell_score, score_hor);
	  cell_score = Max (cell_score, 0);

	  //	  Matrix->grid[i][j] = (Matrix_Cell *) malloc (sizeof(Matrix_Cell));
	  Matrix->grid[i][j]->score =cell_score;
	  if (cell_score == 0)
	    {
	      Matrix->grid[i][j]->prev = NULL;
	      Matrix->grid[i][j]->gap = NOGAP;
	    }
	  else if (cell_score == score_dia)
	    {
	      Matrix->grid[i][j]->prev = Matrix->grid[i-1][j-1];
	      Matrix->grid[i][j]->gap = NOGAP;
	    }
	  else if (cell_score == score_hor)
	    {
	      Matrix->grid[i][j]->prev = Matrix->grid[i][j-1];
	      Matrix->grid[i][j]->gap =H_GAP;
	    }
	  else if (cell_score == score_ver)
	    {
	      Matrix->grid[i][j]->prev = Matrix->grid[i-1][j];
	      Matrix->grid[i][j]->gap =V_GAP;
	    }

	  Matrix->grid[i][j]->row =i;
	  Matrix->grid[i][j]->column =j;

	  if (max_score < cell_score)
	    {
	      temp_node = (Node *) malloc(sizeof(Node));
	      temp_node->row = i;
	      temp_node->column = j;
	      temp_node->prev = NULL;
	      max_cell->size =1;
	      node = max_cell->start;
	      free_Node(node);
	      node = NULL;
	      max_cell->start = temp_node;
	      temp_node = NULL;
	      max_score = cell_score;
	    }
	}
    }
  return max_score;
}

void free_MatrixADT(MatrixADT Matrix)
{
  int i,j;

  for (i =0; i< Matrix->row; i++)
    {
      for (j=0; j< Matrix->column; j++)
	{
	  free(Matrix->grid[i][j]); // (Matrix_Cell *) 
	  Matrix->grid[i][j]= NULL;
	}
      free( Matrix->grid[i]) ; 
      Matrix->grid[i] =NULL;
    }

  free(Matrix);
  Matrix= NULL;
}

void free_Node (Node* node)
{
  if (node == NULL)
    return;
  if ( (node->prev) != NULL) 
    free_Node (node->prev);
  free(node);
  node = NULL;
}

void free_PathADT (PathADT path)
{
  free_Node(path->start);
  path->start = NULL;
  free(path);
  path = NULL;
}













