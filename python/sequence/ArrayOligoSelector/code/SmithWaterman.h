#ifndef _SW
#define _SW

#define DNA_ALPHABET 5
#define UNDEFINED INT_MAX
#define Max(a,b) ( ((a) >= (b)) ? (a) : (b) )  
#define Min(a,b) ( ((a) <= (b)) ? (a) : (b) )  

#define NOGAP 0
#define V_GAP 1
#define H_GAP 2

typedef struct {
  int row;
  int column;
  int **grid;
  int open;
  int extend;
} Score_Matrix;

typedef struct Matrix_Cell{
  int score;
  int gap;
  int row;
  int column;
  struct Matrix_Cell *prev;
} Matrix_Cell;

typedef struct {
  int row;
  int column;
  Matrix_Cell ***grid;
} MatrixCDT;

typedef MatrixCDT *MatrixADT;

typedef struct Node {
  int row;
  int column;
  struct Node *prev;
} Node;

typedef struct PathCDT{
  int size;
  Node *start;
} PathCDT;

typedef PathCDT *PathADT;

enum Matrix_Name {PAM47} ;
enum Grid {ROW, COLUMN};

#endif


