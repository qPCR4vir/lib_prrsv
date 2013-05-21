//core code is from the following : 

/********************************************************************
**
** Copyright (c) 1989 Mark R. Nelson
**
** LZW data compression/expansion demonstration program.
**
** April 13, 1989
**
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BITS 12                   /* Setting the number of bits to 12, 13*/
#define HASHING_SHIFT BITS-8      /* or 14 affects several constants.    */
#define MAX_VALUE (1 << BITS) - 1 /* Note that MS-DOS machines need to   */
#define MAX_CODE MAX_VALUE - 1    /* compile their code in large model if*/
                                  /* 14 bits are selected.               */
#if BITS == 14
#define TABLE_SIZE 18041        /* The string table size needs to be a */
#endif                            /* prime number that is somewhat larger*/
#if BITS == 13                    /* than 2**BITS.                       */
#define TABLE_SIZE 9029
#endif
#if BITS <= 12
#define TABLE_SIZE 5021
#endif

void *malloc();

int *code_value;                  /* This is the code value array        */
unsigned int *prefix_code;        /* This array holds the prefix codes   */
unsigned char *append_character;  /* This array holds the appended chars */
unsigned char decode_stack[4000]; /* This array holds the decoded string */

/********************************************************************
**
** This program gets a file name from the command line.  It compresses the
** file, placing its output in a file named test.lzw.  It then expands
** test.lzw into test.out.  Test.out should then be an exact duplicate of
** the input file.
**
*************************************************************************/
int LZWsize(char* seq, int size) ;
int find_match(int hash_prefix,unsigned int hash_character);
void InputFileToString(FILE* fin, char** seq);

int mainUSEDTOBE(int argc, char *argv[])
{
  FILE* fin;
  int compre_len;
  int size;
  char* sequence;
  int i;
  /*
  **  The three buffers are needed for the compression phase.
  */
  code_value=malloc(TABLE_SIZE*sizeof(unsigned int));
  prefix_code=malloc(TABLE_SIZE*sizeof(unsigned int));
  append_character=malloc(TABLE_SIZE*sizeof(unsigned char));
  if (code_value==NULL || prefix_code==NULL || append_character==NULL)
  {
    printf("Fatal error allocating table space!\n");
    return 0; //error
  }
  /*
  ** Get the file name, open it up, and open up the lzw output file.
  */
  if (argc!=3)
    {
      printf("USAGE: ./LZW sequencefile segmentsize\n");
      return 0;
    }
  else
  {
    /*
    ** Compress the file.
    */
    size = atoi(argv[2]);
    fin = fopen(argv[1],"r");
    if (fin == NULL)
      return 0;
    InputFileToString(fin,&sequence);

    for (i=0; i<strlen(sequence) - size+1 ; i++)
      {
	compre_len = LZWsize(sequence+i,size);
	printf ("%d\n",size - compre_len);
      }
    return 1; //success
  }
}

/*
** This is the compression routine.  The code should be a fairly close
** match to the algorithm accompanying the article.
**
*/

int LZWsize(char* seq, int size) //    (FILE *input,FILE *output)
{
  unsigned int next_code;
  unsigned int character;
  unsigned int string_code;
  unsigned int index;
  int i =0,j=0;

  if (size == 0)
    return 0;
  
  code_value=malloc(TABLE_SIZE*sizeof(unsigned int));
  prefix_code=malloc(TABLE_SIZE*sizeof(unsigned int));
  append_character=malloc(TABLE_SIZE*sizeof(unsigned char));

  next_code=256;              /* Next code is the next available string code*/
  for (i=0;i<TABLE_SIZE;i++)  /* Clear out the string table before starting */
    code_value[i]=-1;

  string_code=seq[0];    /* Get the first code                         */
  i=1;
/*
** This is the main loop where it all happens.  This loop runs util all of
** the input has been exhausted.  Note that it stops adding codes to the
** table after all of the possible codes have been defined.
*/

  while (i<size)
  {
    character = seq[i];
    ++i;

    index=find_match(string_code,character);/* See if the string is in */
    if (code_value[index] != -1)            /* the table.  If it is,   */
      string_code=code_value[index];        /* get the code value.  If */
    else                                    /* the string is not in the*/
    {                                       /* table, try to add it.   */
      if (next_code <= MAX_CODE)
      {
        code_value[index]=next_code++;
        prefix_code[index]=string_code;
        append_character[index]=character;
      }
      // output_code(output,string_code);  /* When a string is found  */
      string_code=character;            /* that is not in the table*/
      ++j;
    }                                   /* I output the last string*/
  }                                     /* after adding the new one*/
/*
** End of the main loop.
*/

  ++j;

  free(code_value);
  free(prefix_code);
  free(append_character);

  return j;
}

/*
** This is the hashing routine.  It tries to find a match for the prefix+char
** string in the string table.  If it finds it, the index is returned.  If
** the string is not found, the first available index in the string table is
** returned instead.
*/

int find_match(int hash_prefix,unsigned int hash_character)
{
  int index;
  int offset;

  index = (hash_character << ( HASHING_SHIFT )) ^ hash_prefix ;
  if (index == 0)
    offset = 1;
  else
    offset = TABLE_SIZE - index;
  while (1)
  {
    if (code_value[index] == -1)
      return(index);
    if (prefix_code[index] == hash_prefix && 
        append_character[index] == hash_character)
      return(index);
    index -= offset;
    if (index < 0)
      index += TABLE_SIZE;
  }
}

/*
**  This is the expansion routine.  It takes an LZW format file, and expands
**  it to an output file.  The code here should be a fairly close match to
**  the algorithm in the accompanying article.
*/

/*
** This routine simply decodes a string from the string table, storing
** it in a buffer.  The buffer can then be output in reverse order by
** the expansion program.
*/

char *decode_string(unsigned char *buffer,unsigned int code)
{
int i;

  i=0;
  while (code > 255)
  {
    *buffer++ = append_character[code];
    code=prefix_code[code];
    if (i++>=4094)
    {
      printf("Fatal error during code expansion.\n");
      exit(0);
    }
  }
  *buffer=code;
  return(buffer);
}
