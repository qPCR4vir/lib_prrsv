#include <stdlib.h>
#include <string.h>
#include <stdio.h>

void InputFileToString(FILE* fin, char** seqPtr)
{
  int seqSize = 1024;
  char c;
  int count =0;
  char* newPtr;

  *seqPtr = (char*) malloc(sizeof(char)*1024);
  memset(*seqPtr,'\0',1024);

  do 
    {
      c = getc(fin);
      if (c== EOF)
	break;
      (*seqPtr)[count] =c;
      count++;
      if (count >= seqSize -2)
	{
	  seqSize += 1024;
	  (*seqPtr)[count] ='\0';
	  newPtr = realloc(*seqPtr,seqSize);
	  *seqPtr = newPtr;
	}
    }
  while (!feof(fin));
  (*seqPtr)[count] ='\0';
  fclose(fin);
}

