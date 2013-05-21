#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

char Complement_Letter_DRNA ( char c)
{
  switch (toupper(c))
    {
    case 'A':
      return 'T';
    case 'T':
      return 'A';
    case 'G':
      return 'C';
    case 'C':
      return 'G';
    case 'U':
      return 'A';
    default:
      return 'N';
    }
}

void Complement_Seq_RDNA (char *seq, char *comp_seq)
{
  int i;
  int l = strlen(seq);

  for (i =0 ; i < l; i++)
    {
      comp_seq[i] = Complement_Letter_DRNA (seq[i]);
    }
}
	    










