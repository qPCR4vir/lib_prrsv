//=============================================================================
// Module:        galign.cpp
// Project:       Cubic Project - Calculation of melting temperature and free
//                energy of two DNA strands 
// Type:          implementation - Thermodynamic Alignment.
// Language:      c++
// Compiler:      microsoft visual c++ 6.0, unix/linux gcc
// System/OS:     Windows 32, Sun solaris, Linux, other unix systems (untested)
// Database:      none
// Description:   class GAlign - Thermodynamic Alignment Algorithm
// Author:        leber
// Date:          01/2002 - 02/2002
// Copyright:     (c) L. Kaderali & M. Leber, 01/2002 - 02/2002
//
// Revision History
// $ 00sep04 LK : created 
// $ 00dec30 LK : changed to do local alignment of probe against
//                one entire sequence 
// $ 01feb07 LK : optimized!
// $ 01aug06 LK : corrected, included salt and concentration input;
// $              made true local alignment (problem with initial mismatches!) 
// #$
//=============================================================================

#include <math.h>
#include <assert.h>
#include <iomanip>
#include <stdio.h>
#ifdef _DEBUG
#include <fstream>
#endif
#include "galign.h"

#ifdef _pack
  #pragma pack(1)
  #define dH(a,b,c) (*(dH+((b)*(maxseq1len+1)*3+(a)*3+(c))))
  #define dS(a,b,c) (*(dS+((b)*(maxseq1len+1)*3+(a)*3+(c))))
#else
  #define dH(a,b,c) (*(dH+((b)*(maxseq1len+1)*3+(a)*3+(c))))
  #define dS(a,b,c) (*(dS+((b)*(maxseq1len+1)*3+(a)*3+(c))))
#endif



///////////////////////////////////////////////////////////////////////////////
// Construction/Destruction
///////////////////////////////////////////////////////////////////////////////

//------------------------ Constructor ----------------------------------------
// Initialize Thermalign Object. This will allocate memory for the dynamic
// programming matrices for both entropy and enthalpy of the present
// alignment. The max length of the two sequences to align has to be passed 
// to the routine. Also, an initialized CNNParams Object must be passed to 
// CThermAlign. Object will be initialized only once, sequences can be
// passed to it later.

//int s1lmax,s2lmax;  // global requiered for DEBUGGING

GAlign::GAlign(int ALength, int BLength, PNNParams myNNParams)
{
    // Create dynamic programming matrices for entropy and enthalpy.
    // The first two dimensions are for the two strings respectively, the
    // third dimensions differentiates the alignment according to its end:
    // type 0 Alignments will be those that align x(i+1) with y(j+1),
    // type 1 Alignments align x(i+1) with a gap
    // type 2 Alignments align y(i+1) with a gap
    // All three of these types must separately be considered to select the
    // optimum alignment in the next iteration.
    // The alignment algorithm proceeds analogous to the well known
    // Smith-Waterman algorithm for global alignments.

    dH = new float[(ALength+1)*(BLength+1)*3];
    dS = new float[(ALength+1)*(BLength+1)*3];
//    s1lmax=ALength;
//    s2lmax=BLength;

	if ((dH==NULL)||(dS==NULL))
	{
		printf("Fatal Error: Out of memory!\n");
		exit(-1);
	}

    GNNParams = myNNParams;
	seq1len = ALength;
	maxseq1len = seq1len;
	seq2len = BLength;

    // set dH / dS entries all to zero!
    #ifdef _pack
	    // this is FAST, but works only if #pragma pack(1) has been set...
	    memset(dH,0,(ALength+1)*(BLength+1)*3*sizeof(float));
	    memset(dS,0,(ALength+1)*(BLength+1)*3*sizeof(float));
    #else
		// the following is slower, but works also when padding bytes have 
		// been inserted in the array for alignment by the compiler...
		// time is not really critical here, as this initialization is done
		// only once...
		for (int a=0;a<=ALength;a++)
			for (int b=0;b<=BLength;b++)
			{
				dH(a,b,0)=0.0f;
				dH(a,b,1)=0.0f;
				dH(a,b,2)=0.0f;
				dS(a,b,0)=0.0f;
				dS(a,b,1)=0.0f;
				dS(a,b,2)=0.0f;
			}
    #endif
    InitBorder();  // also need to do only once...
//    forbidden_entropy = (-(GNNParams->R)*(float)log((GNNParams->Ct)/4.0f));
    return;
}

//------------------------ Destructor -----------------------------------------
// Free memory

GAlign::~GAlign()
{
	delete[] dH;
	delete[] dS;
    return;
}



///////////////////////////////////////////////////////////////////////////////
//  FUNCTION:   void GAlign::InitStrings(char* s1, 
//                                            char* s2,
//                                            int s1len,
//                                            int s2len);
//
//  PURPOSE:    Initialize strings for alignment table!
//
//  PARAMETERS:
//      s1:    Target sequence
//      s2:    Probe sequence
//      s1len: Target sequence length
//      s2len: Probe sequence length
//
//  RETURN VALUE:
//      void
//
//  REVISION HISTORY
//  $ 01jan03 LK : created
//  $ 01feb07 LK : removed maximum local stuff
//  #$

void GAlign::InitStrings(char* s1,char* s2,int s1len ,int s2len,int minxpos)
{
	seq1=s1; seq1len=s1len;
	seq2=s2; seq2len=s2len;
	GNNParams -> UpdateParams(s1,s2);
/* lk01feb07: removed maxloc-stuff
	// if maximum local is not in reused subtable, need to reset!
	if ((maxlocj>=minxpos)||(maxloci!=seq1len))
	{
		maxloci=0; maxlocj=0; maxlocg=0;
	}
*/
//    InitBorder();   //lk01feb07: removed, as it suffices if this is done once in constructor
	return;
}


///////////////////////////////////////////////////////////////////////////////
//  FUNCTION:   void GAlign::InitBorder();
//
//  PURPOSE:    Initialize the lower and left border of the dynamic
//              programming table. 
//
//  PARAMETERS:
//      none
//
//  RETURN VALUE:
//      void
//
//  REVISION HISTORY
//  $ 00sep06 LK : created
//  $ 00dec30 LK : complete revision - simply init table now, no
//                 more saving etc...
//  $ 01feb07 LK : removed local alignment stuff and optimized
//  #$

void GAlign::InitBorder()
{
	// Initialization for new alignment!
	maxlocg = -99999;   // maximum melting temperature found (kelvin)
	maxloci  = 0;
	maxlocj  = 0;
	maxloct  = 0;

	for (int i=0; i<=seq1len; i++)
        {
		dH(i,0,0) = 0.0f;	dH(i,0,1) = 0.0f;	dH(i,0,2) = 0.0f;
		dS(i,0,0) = GNNParams->GetInitialEntropy();	
                dS(i,0,1) = GNNParams->GetInitialEntropy();	
                dS(i,0,2) = GNNParams->GetInitialEntropy();     
        }     
        for  (int j=1; j<=seq2len; j++)     {
		dH(0,j,0) = 0.0f;	dH(0,j,1) = 0.0f;	dH(0,j,2) = 0.0f;
		dS(0,j,0) = GNNParams->GetInitialEntropy(); 
                dS(0,j,1) = GNNParams->GetInitialEntropy();
                dS(0,j,2) = GNNParams->GetInitialEntropy();  
        }

    return;
}


///////////////////////////////////////////////////////////////////////////////
//  FUNCTION:   GAlign::CalculateTable();
//
//  PURPOSE:    Compute the entire dynamic programming table cell by cell.
//              Requieres that InitBorder has been called.
//              Computations will be done in the usual rowwise manner.
//
//  PARAMETERS:
//      int startRow:  Row that computation shall commence at!
//
//  RETURN VALUE:
//      void - Data will be written directly to dH and dS in this object 
//
//  REVISION HISTORY
//  $              00sep06 : created LK
//                 00dec31 : modified to begin at a specified row!
//  #$

void GAlign::CalculateTable(int startRow)
{
//	#define forbidden_enthalpy 1000000000000000000.0f
	register int i;
	register int iminusone;
	register int j;
	register int jminusone;
	// Determine previous character!
	register char prevchari;
	register char prevcharj;
	register char seq1char;
	register char seq2char;
	register float entropy1;
	register float entropy2;
	register float entropy3;
	register float enthalpy1;
	register float enthalpy2;
	register float enthalpy3;
	register float g1;
	register float g2;
	register float g3;
	float maxg;
//    ofstream dboutstream("debug.out");
//    assert(startRow>0);
    for (j=startRow; j<=seq2len; j++)
	{
		if (j>1)
			prevcharj = seq2[j-2];
		else
			prevcharj = 0;
		jminusone = j - 1;
		seq2char = seq2[jminusone];
        for (i=1; i<=seq1len; i++)
        {
		   // local variables
			if (i>1)
				prevchari = seq1[i-2];
			else
				prevchari = 0;
			iminusone = i - 1;
			seq1char=seq1[iminusone];

			// --------------- Upper cell: Alignment of seq1[i] with seg2[j]
			entropy1 = dS(iminusone,jminusone,0);
			enthalpy1 = dH(iminusone,jminusone,0) ;
			entropy2 = dS(iminusone,jminusone,1);
			enthalpy2 = dH(iminusone,jminusone,1);
			entropy3 = dS(iminusone,jminusone,2);
			enthalpy3 = dH(iminusone,jminusone,2);
			if (enthalpy1<forbidden_enthalpy)
			{
			     entropy1 += GNNParams->GetEntropy(prevchari,seq1char, prevcharj, seq2char);
			     enthalpy1 +=GNNParams->GetEnthalpy(prevchari, seq1char, prevcharj, seq2char);

		        } 			
                        if (enthalpy2<forbidden_enthalpy)
			{
			     enthalpy2 +=GNNParams->GetEnthalpy(prevchari,seq1char, 0, seq2char); // 0 is gap!
                             entropy2 += GNNParams->GetEntropy(prevchari,seq1char, 0, seq2char);
 			} 			
                        if(enthalpy3<forbidden_enthalpy) 			
                        {
		             enthalpy3 += GNNParams->GetEnthalpy(0, seq1char,prevcharj, seq2char);
                             entropy3 +=GNNParams->GetEntropy(0, seq1char, prevcharj, seq2char);
 	
                  	} 			
                        // choose optimum combination
			g1 = GNNParams->CalcG(entropy1, enthalpy1);
			g2 = GNNParams->CalcG(entropy2, enthalpy2);
			g3 = GNNParams->CalcG(entropy3, enthalpy3);
			if ((g1>=g2)&&(g1>=g3))
			{
				dS(i,j,0) = entropy1;
				dH(i,j,0) = enthalpy1;
				maxg=g1;
			}
			else if (g2>=g3)
			{
				dS(i,j,0) = entropy2;
				dH(i,j,0) = enthalpy2;
				maxg=g2;
			}
			else
			{
				dS(i,j,0) = entropy3;
				dH(i,j,0) = enthalpy3;
				maxg=g3;
			}
		    // set to zero if alignment so far is mismatches only (tm will be zero -> no need to reset maxg)!
                        if(dH(i,j,0)==0)
                            dS(i,j,0)=GNNParams->GetInitialEntropy();
			// check if local maximum found!
			if (((i==seq1len)||(j==seq2len))&&(maxg>=maxlocg))
			{
				maxlocg=maxg;
				maxloci=i;
				maxlocj=j;
				maxloct=0;
			}
		    // set to zero if alignment so far is mismatches only!
			if (dH(i,j,0)==0)
			    dS(i,j,0)=GNNParams->GetInitialEntropy();

			// --------------- Middle cell: Alignment of seq1[i] with '-'
			// following changes lk 01jan08
			// we do not allow $- in the beginning, therefore set to forbidden if
			// j=1
			if (j==1)
			{
				enthalpy1=forbidden_enthalpy; enthalpy2=forbidden_enthalpy;
				entropy1=forbidden_entropy; entropy2=forbidden_entropy;
			}
			// also, disallow -$ in the end, therefore forbid if j=seq2len-1
			else if (j==seq2len-1)
			{
				enthalpy1=forbidden_enthalpy; enthalpy2=forbidden_enthalpy;
				entropy1=forbidden_entropy; entropy2=forbidden_entropy;
			}
			else
			{
				// end lk01jan08
				// -- entropy
				entropy1 = dS(iminusone,j,0);
				entropy2 = dS(iminusone,j,1);
				enthalpy1 = dH(iminusone,j,0);
				enthalpy2 = dH(iminusone,j,1);
				if (enthalpy1<forbidden_enthalpy)
				{
				    entropy1+=GNNParams->GetEntropy(prevchari, seq1char, seq2char, 0);
 	                            enthalpy1 +=GNNParams->GetEnthalpy(prevchari,seq1char, seq2char, 0);
                                }  	 
                                if (enthalpy2<forbidden_enthalpy)
                                {
	 	 		    entropy2 +=GNNParams->GetEntropy(prevchari, seq1char, 0, 0); 
                                    enthalpy2 +=GNNParams->GetEnthalpy(prevchari, seq1char, 0, 0); 		
	 			}
                        }	 	
                        // -- choose optimum combination 
                        g1 = GNNParams->CalcG(entropy1, enthalpy1);
                        g2 = GNNParams->CalcG(entropy2, enthalpy2);  
                        if (g1>=g2)
                        { 	
                                    dS(i,j,1)=entropy1; 	
                                    dH(i,j,1)=enthalpy1; 	
                                    maxg=g1;
                        }
			else
			{
				dS(i,j,1)=entropy2;
				dH(i,j,1)=enthalpy2;
				maxg=g2;
			}
			// check if local maximum found!
			if (((i==seq1len)||(j==seq2len))&&(maxg>=maxlocg))
			{
				maxlocg=maxg;
				maxloci=i;
				maxlocj=j;
				maxloct=1;
			}
		    // set to zero if alignment so far is mismatches only!
			if (dH(i,j,1)==0)
			    dS(i,j,1)=GNNParams->GetInitialEntropy();
		
			// --------------- Lower cell: Alignment of '-' with seq2[j]
			// following changes lk 01jan08
			// we do not allow $- in the beginning, therefore set to forbidden if
			// i=1
			if (i==1)
			{
				enthalpy1=forbidden_enthalpy; enthalpy2=forbidden_enthalpy;
				entropy1=forbidden_entropy; entropy2=forbidden_entropy;
			}
			// also, disallow -$ in the end, therefore forbid if i=seq1len-1
			else if (i==seq1len-1)
			{
				enthalpy1=forbidden_enthalpy; enthalpy2=forbidden_enthalpy;
				entropy1=forbidden_entropy; entropy2=forbidden_entropy;
			}
			else
			{
			// end lk01jan08
				entropy1 = dS(i,jminusone,0);
				entropy2 = dS(i,jminusone,2);
				enthalpy1 = dH(i,jminusone,0);
				enthalpy2 = dH(i,jminusone,2);
				if (enthalpy1<forbidden_enthalpy)
				{
					entropy1+=GNNParams->GetEntropy(seq1char, 0, prevcharj, seq2char);
 				        enthalpy1 +=GNNParams->GetEnthalpy(seq1char, 0, prevcharj, seq2char); 
                                } 		 	
                                if(enthalpy2<forbidden_enthalpy)
 				{ 	
                                        entropy2 += GNNParams->GetEntropy(0,0, prevcharj, seq2char); 
                                        enthalpy2 +=GNNParams->GetEnthalpy(0, 0, prevcharj, seq2char);
                                } 	 	
                        } 			
                        // -- choose optimum combination 	
                        g1 = GNNParams->CalcG(entropy1, enthalpy1); 	
                        g2 = GNNParams->CalcG(entropy2, enthalpy2); 	
	                if (g1>=g2)
                        {
				dS(i,j,2)=entropy1;
				dH(i,j,2)=enthalpy1;
				maxg=g1;
			}
			else
			{
				dS(i,j,2)=entropy2;
				dH(i,j,2)=enthalpy2;
				maxg=g2;
			}
			// check if local maximum found!
			if (((i==seq1len)||(j==seq2len))&&(maxg>=maxlocg))
			{
				maxlocg=maxg;
				maxloci=i;
				maxlocj=j;
				maxloct=2;
			}
		    // set to zero if alignment so far is mismatches only!
			if (dH(i,j,2)==0)
			    dS(i,j,2)=GNNParams->GetInitialEntropy();
			
		
		/*
            dboutstream<<"After cell ("<<i<<","<<j<<"):"<<endl;
            PrintDPTable(dboutstream);
            dboutstream<<endl<<flush;
		*/
        }
	}
    //printEnthalpyTable(0);
    //printEntropyTable(0);
    return;
}

void GAlign::printEnthalpyTable(int level)
{
	for(int i=0;i<=seq1len;i++) 
	{
		for(int j=0;j<=seq2len;j++)
		{
			cout << dH(i,j,level) << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void GAlign::printEntropyTable(int level)
{
	for(int i=0;i<=seq1len;i++) 
	{
		for(int j=0;j<=seq2len;j++)
		{
			cout << dS(i,j,level) << " ";
		}
		cout << endl;
	}
	cout << endl;
}



#ifdef _output_alignment
///////////////////////////////////////////////////////////////////////////////
//  FUNCTION:   void GAlign::OutputLocalAlignment(outputStream);
//
//  PURPOSE:    Output optimum alignment (Local). Requieres that 
//              CalculateTable() has finished.
//
//  PARAMETERS:
//      none
//
//  RETURN VALUE:
//      float
//
//  REVISION HISTORY
//  $              00sep06 : created LK & ML
//  $              00nov16 : changed to print to outputStream
//  $              03jan26 :
//  $              03jan26 : changed return value to float
//  #$

void GAlign::OutputLocalAlignment(ostream &outputStream)
{
	s1align.freeze(0);
	s2align.freeze(0);
	aligntype.freeze(0);
	s1align.seekp(0);
	s2align.seekp(0);
	aligntype.seekp(0);

	//outputStream<<"Maximum Local Alignment: " << endl;

	s1aptr=&s1align;
	s2aptr=&s2align;
	atypptr=&aligntype;
	int maxloct;
        float TempK;

	float g0, g1, g2;
	g0 = GNNParams->CalcG(dS(maxloci,maxlocj,0),dH(maxloci,maxlocj,0));
	g1 = GNNParams->CalcG(dS(maxloci,maxlocj,1),dH(maxloci,maxlocj,1));
	g2 = GNNParams->CalcG(dS(maxloci,maxlocj,2),dH(maxloci,maxlocj,2));

	if ((g0>=g1)&&(g0>=g2))
	{
   	    //outputStream<<"Free energy = "<< g0 <<endl;
		maxloct=0;
	}
	else if ((g1>=g0)&&(g1>=g2))
	{
   	    //outputStream<<"Free energy = "<< g1 <<endl;
		maxloct=1;
	}
	else
	{
   	    //outputStream<<"Free energy = "<< g2 <<endl;
		maxloct=2;
	}

	bool good = OutputAlignment(outputStream, maxloci, maxlocj, maxloct, true);
	s1align<<'\0';
	s2align<<'\0';
	aligntype<<'\0';
	if (good)
        {
		outputStream<<s1align.str()<<endl<<s2align.str()<<endl <<aligntype.str()<<endl<<flush;
        }
	else
		outputStream<<"Alignment Error!"<<endl<<flush;

    return;
}
#endif


#ifdef _output_alignment
///////////////////////////////////////////////////////////////////////////////
//  FUNCTION:   void GAlign::OutputAlignment(int i, int j, int t, bool local=false); 
//  PURPOSE:    Output optimum alignment of prefixes seq1[1..i] and seq2[1..j].
//              Requieres that CalculateTable() has finished.
//              Procedes recursively by calling self for shorter sequences 
//              and adding the final character.
//
//  PARAMETERS:
//      int i   - end of prefix of sequence 1 to output
//      int j   - end of prefix of sequence 2 to output
//      int t   - type of alignment:
//                 0 : (i,j) aligned
//                 1 : (i,-) aligned
//                 2 : (-,j) aligned
//      bool local - set to true if this is a local alignment 
//
//  RETURN VALUE:
//      bool    - true if backtracking path is valid, false otherwise.
//                (Prints directly to screen)
//
//  REVISION HISTORY
//  $              00sep06 : created LK
//  $              00nov16 : changed to output to file
//  #$

bool GAlign::OutputAlignment(ostream &outputStream, int i, int j, int t, bool local) {
    // have to backtrace through table to find optimum alignment!
    // do this recursively by calling optimum alignment for prefixes of the
    // present alignment...
    #define matchchar "M"
    #define mismatchchar "."
    #define insertchar "."
    #define deletechar "."
	bool good = false;
    char prevchari = '*';
    char prevcharj = '*';
    if (i>1)
        prevchari = seq1[i-2];
    if (j>1)
        prevcharj = seq2[j-2];

	// check if done
	if ((i<=0) && (j<=0))
		return true;  // done!

	// if at border: initial gaps incur no cost! -> done!
	if (i==0)
	{
	    if (local==false)
		    OutputAlignment(outputStream,i,j-1,2);
		(*s1aptr)<<"-";
		(*s2aptr)<<seq2[j-1];
		(*atypptr)<<deletechar;
		return true;
	}
	if (j==0)
	{
	    if (local==false)
		  OutputAlignment(outputStream,i-1,j,1);
		(*s1aptr)<<seq1[i-1];
		(*s2aptr)<<"-";
		(*atypptr)<<insertchar;
		return true;
	}

	// which type is the best alignment?
	float h0, h1, h2, s0, s1, s2;
	if (t==0)   // align (i,j). What was previous alignment?
	{
		h0= dH(i-1,j-1,0) 
			+GNNParams->GetEnthalpy(prevchari,seq1[i-1],prevcharj,seq2[j-1]);
		h1= dH(i-1,j-1,1) 
			+ GNNParams->GetEnthalpy(prevchari,seq1[i-1],'-',seq2[j-1]);
		h2= dH(i-1,j-1,2) 
			+ GNNParams->GetEnthalpy('-',seq1[i-1],prevcharj,seq2[j-1]);
		s0= dS(i-1,j-1,0) 
			+ GNNParams->GetEntropy(prevchari,seq1[i-1],prevcharj,seq2[j-1]);
		s1= dS(i-1,j-1,1) 
			+ GNNParams->GetEntropy(prevchari,seq1[i-1],'-',seq2[j-1]);
		s2= dS(i-1,j-1,2) 
			+ GNNParams->GetEntropy('-',seq1[i-1],prevcharj,seq2[j-1]);
		if ((h0==dH(i,j,0))&&(s0==dS(i,j,0)))
			good=OutputAlignment(outputStream,i-1,j-1,0,local);
		if ((h1==dH(i,j,0))&&(s1==dS(i,j,0))&&(!good))
			good=OutputAlignment(outputStream,i-1,j-1,1,local);
		if ((h2==dH(i,j,0))&&(s2==dS(i,j,0))&&(!good))
			good=OutputAlignment(outputStream,i-1,j-1,2,local);
		if (good)
		{
			(*s1aptr)<<seq1[i-1];
			(*s2aptr)<<seq2[j-1];
			if (seq1[i-1]=='A')
			{
				if (seq2[j-1]=='T')
					(*atypptr)<<matchchar;
				else
					(*atypptr)<<mismatchchar;
			}
			if (seq1[i-1]=='C')
			{
				if (seq2[j-1]=='G')
					(*atypptr)<<matchchar;
				else
					(*atypptr)<<mismatchchar;
			}
			if (seq1[i-1]=='G')
			{
				if (seq2[j-1]=='C')
					(*atypptr)<<matchchar;
				else
					(*atypptr)<<mismatchchar;
			}
			if (seq1[i-1]=='T')
			{
				if (seq2[j-1]=='A')
					(*atypptr)<<matchchar;
				else
					(*atypptr)<<mismatchchar;
			}
			if (seq1[i-1]=='$')
				(*atypptr)<<mismatchchar;
		}
	}
	else if (t==1)  // align (i,-)
	{
		h0= dH(i-1,j,0) 
			+ GNNParams->GetEnthalpy(prevchari,seq1[i-1],seq2[j-1],'-');
		h1= dH(i-1,j,1) 
			+ GNNParams->GetEnthalpy(prevchari,seq1[i-1],'-','-');
		s0= dS(i-1,j,0) 
			+ GNNParams->GetEntropy(prevchari,seq1[i-1],seq2[j-1],'-');
		s1= dS(i-1,j,1) 
			+ GNNParams->GetEntropy(prevchari,seq1[i-1],'-','-');
		if ((h0==dH(i,j,1))&&(s0==dS(i,j,1)))
			good= OutputAlignment(outputStream,i-1,j,0,local);
		if ((h1==dH(i,j,1))&&(s1==dS(i,j,1))&&(!good))
			good= OutputAlignment(outputStream,i-1,j,1,local);
        if (good)
		{
			(*s1aptr)<<seq1[i-1];
			(*s2aptr)<<"-";
			(*atypptr)<<insertchar;
		}
	}
	else if (t==2)
	{
		h0= dH(i,j-1,0) 
			+ GNNParams->GetEnthalpy(seq1[i-1],'-',prevcharj,seq2[j-1]);
		h2= dH(i,j-1,2) 
			+ GNNParams->GetEnthalpy('-','-',prevcharj,seq2[j-1]);
		s0= dS(i,j-1,0) 
			+ GNNParams->GetEntropy(seq1[i-1],'-',prevcharj,seq2[j-1]);
		s2= dS(i,j-1,2) 
			+ GNNParams->GetEntropy('-','-',prevcharj,seq2[j-1]);
		if ((h0==dH(i,j,2))&&(s0==dS(i,j,2)))
			good = OutputAlignment(outputStream,i,j-1,0,local);
		if ((h2==dH(i,j,2))&&(s2==dS(i,j,2))&&(!good))
			good = OutputAlignment(outputStream,i,j-1,2,local);
		if (good)
		{
			(*s1aptr)<<"-";
			(*s2aptr)<<seq2[j-1];
			(*atypptr)<<deletechar;
		}
	}
	if (local)
	{
	    if (!good)
	    {
		    if ((i>0)&&(j>0))
			{
			    (*s1aptr)<<seq1[i-1];
				(*s2aptr)<<seq2[j-1];
				(*atypptr)<<",";
			}
	    }
		return true;
	}
	else
	  return good;
}
#endif


///////////////////////////////////////////////////////////////////////////////
//  FUNCTION:   float GetEntropy(int i, int j);
//
//  PURPOSE:    Return Entropy for given cell
//
//  PARAMETERS:
//      i, j    - Cell coordinates
//
//  RETURN VALUE:
//      float  - Entropy (dS) 
//
//  REVISION HISTORY
//  $              00sep23 : created LK
//  #$

/* inline */ float GAlign::GetEntropy(int i, int j)
{
	float entropy;

	// first, need to determine optimum values for dG and dH
	float tm0, tm1, tm2;
	tm0 = GNNParams->CalcTM(dS(i,j,0),dH(i,j,0));
	tm1 = GNNParams->CalcTM(dS(i,j,1),dH(i,j,1));
	tm2 = GNNParams->CalcTM(dS(i,j,2),dH(i,j,2));
	if ((tm0>=tm1)&&(tm0>=tm2))
		entropy = dS(i,j,0);
	else if ((tm1>=tm0)&&(tm1>=tm2))
		entropy = dS(i,j,1);
	else
		entropy = dS(i,j,2);

	return entropy;
}


///////////////////////////////////////////////////////////////////////////////
//  FUNCTION:   float GetEnthalpy(int i, int j);
//
//  PURPOSE:    Return Enthalpy for given cell
//
//  PARAMETERS:
//      i, j    - Cell coordinates
//
//  RETURN VALUE:
//      float  - Enthalpy (dS) 
//
//  REVISION HISTORY
//  $              00sep23 : created LK
//  #$

/* inline */ float GAlign::GetEnthalpy(int i, int j)
{
	float enthalpy;

	// first, need to determine optimum values for dG and dH
	float tm0, tm1, tm2;
	tm0 = GNNParams->CalcTM(dS(i,j,0),dH(i,j,0));
	tm1 = GNNParams->CalcTM(dS(i,j,1),dH(i,j,1));
	tm2 = GNNParams->CalcTM(dS(i,j,2),dH(i,j,2));
	if ((tm0>=tm1)&&(tm0>=tm2))
		enthalpy = dH(i,j,0);
	else if ((tm1>=tm0)&&(tm1>=tm2))
		enthalpy = dH(i,j,1);
	else
		enthalpy = dH(i,j,2);

	return enthalpy;
}



///////////////////////////////////////////////////////////////////////////////
//  FUNCTION:   float GetFreeEnergyK(int i, int j, float t);
//
//  PURPOSE:    Return free Energy for given cell at given temperature
//
//  PARAMETERS:
//      i, j    - Cell coordinates
//      t       - Temperature in Kelvin
//
//  RETURN VALUE:
//      float  - free Energy value (dG) at given temperature (Kelvin)
//
//  REVISION HISTORY
//  $              00sep06 : created LK
//  #$

float GAlign::GetFreeEnergyK(int i, int j, float t)
{
	// dG = dH - T * dS 

	float entropy;
	float enthalpy;

	// first, need to determine optimum values for dG and dH
	entropy = -GetEntropy(i,j);
	enthalpy = -GetEnthalpy(i,j);

	return enthalpy - t * entropy;
}


///////////////////////////////////////////////////////////////////////////////
//  FUNCTION:   float GetFreeEnergyC(int i, int j, int t);
//
//  PURPOSE:    Return free Energy for given cell at given temperature
//
//  PARAMETERS:
//      i, j    - Cell coordinates
//      t       - Temperature in Celsius
//
//  RETURN VALUE:
//      float  - free Energy value (dG) at melting temperature
//
//  REVISION HISTORY
//  $              00sep06 : created LK
//  #$

float GAlign::GetFreeEnergyC(int i, int j, float t)
{
	// dG = dH - T * dS

	float entropy;
	float enthalpy;

	// first, need to determine optimum values for dG and dH
	entropy = GetEntropy(i,j);
	enthalpy = GetEnthalpy(i,j);

	return enthalpy - (t + 273.15f) * entropy;
}

///////////////////////////////////////////////////////////////////////////////
//  FUNCTION:   float GAlign::GetMeltingTempC(int i, int j)
//
//  PURPOSE:    Calculate TM for given cell. Return in Celsius.
//
//  PARAMETERS:
//      i, j    - Cell coordinates
//
//  RETURN VALUE:
//      float  - Melting Temperature TM
//
//  REVISION HISTORY
//  $              00sep06 : created LK
//  #$

float GAlign::GetMeltingTempC(int i, int j)
{
	return GetMeltingTempK(i,j) - 273.15f;
}


///////////////////////////////////////////////////////////////////////////////
//  FUNCTION:   float GAlign::GetMeltingTempK(int i, int j)
//
//  PURPOSE:    Calculate TM for given cell. Return in Kelvin.
//
//  PARAMETERS:
//      i, j    - Cell coordinates
//
//  RETURN VALUE:
//      float  - Melting Temperature TM
//
//  REVISION HISTORY
//  $              00sep06 : created LK
//  #$

float GAlign::GetMeltingTempK(int i, int j)
{
	float tm0, tm1, tm2;
	tm0 = GNNParams->CalcTM(dS(i,j,0),dH(i,j,0));
	tm1 = GNNParams->CalcTM(dS(i,j,1),dH(i,j,1));
	tm2 = GNNParams->CalcTM(dS(i,j,2),dH(i,j,2));
	float maxg;
	if ((tm0>=tm1)&&(tm0>=tm2))
		maxg=tm0;
	else if ((tm1>=tm0)&&(tm1>=tm2))
		maxg=tm1;
	else
		maxg=tm2;
	if (maxg<0)
		maxg=0;
	return maxg;
}


