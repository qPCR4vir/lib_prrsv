//=============================================================================
// Module:        galign.h
// Project:       Cubic Project - Calculation of melting temperature and free
//                energy of two DNA strands 
// Type:          header file - Thermodynamic Alignment.
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
// $ 01feb07 LK : optimized
// #$
//=============================================================================

#if !defined(AFX_GAlign_H__1B9227F7_82AB_11D4_9FFF_000000000000__INCLUDED_)
#define AFX_GAlign_H__1B9227F7_82AB_11D4_9FFF_000000000000__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <iostream>
#include <strstream>
//#include <ostream>
#include "nnparams.h"                         // Nearest Neighbor Parameters
//#include "../gsfxtree/gsfxtree.h"             // Suffix Tree Stuff

using namespace std;

#ifdef _pack
#pragma pack(1)
#endif

//-----------------------------------------------------------------------------
// class GAlign

typedef class GAlign* GGAlign;

class GAlign  
{
public:
	GAlign(int, int, PNNParams);
	virtual ~GAlign();
	void InitStrings(char*,char*,int,int,int=0);
	void InitBorder();
	void CalculateTable(int=1);
	#ifdef _output_alignment
		bool OutputAlignment(ostream &outputStream, int, int, int, bool local=false);
		void OutputLocalAlignment(ostream &outputStream);
	#endif
	float GetEntropy(int, int);
	float GetEnthalpy(int, int);
	float GetFreeEnergyK(int, int, float);
	float GetFreeEnergyC(int, int, float);
	float GetMeltingTempC(int, int);
	float GetMeltingTempK(int, int);

	void printEnthalpyTable(int level);
	void printEntropyTable(int level);

// lk01feb07: removed maxlocstuff as not requiered by thermtreealign...
	float maxlocg;      // maximum local dG value found
	int maxloci;         // i position thereof
	int maxlocj;         // j position thereof
	int maxloct;         // and type of maximum alignment!
    int targetNumber;    // id number of target sequence
//    PGSfxLeaf probeNode; // identifier for probe sequence
private:
	char *seq1;         // Sequence 1
	char *seq2;         // Sequence 2
	int seq1len;        // Length of Sequence 1
	int seq2len;        // Length of Sequence 2
	int maxseq1len;     // length of longest target...
	float *dH;         // Dynamic Programming Table for Entropy
	float *dS;         // Dynamic Programming Table for Enthalpy
	PNNParams GNNParams; // Nearest Neighbor parameters
//	float forbidden_entropy;
#ifdef _output_alignment
	ostrstream *s1aptr, *s2aptr, *atypptr;
// Used to buffer aligned sequences (on output)
	ostrstream s1align;
	ostrstream s2align;
	ostrstream aligntype;   // insert, deletion, match, unmatch (for output)
	// same for local alignment:
	/*	removed lk00jan08: use global instead!
	ostrstream ls1align;   
	ostrstream ls2align;
	ostrstream laligntype;  // insert, deletion, match, unmatch (for output)
	*/
#endif
};

#endif //!defined(AFX_GAlign_H__1B9227F7_82AB_11D4_9FFF_000000000000__INCLUDED_)

