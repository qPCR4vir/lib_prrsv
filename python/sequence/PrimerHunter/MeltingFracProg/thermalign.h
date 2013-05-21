//=============================================================================
// Module:        thermalign.h
// Project:       Diploma Thesis - Probe Selection for DNA Microarrays
// Type:          header file - Thermodynamic Alignment.
// Language:      c++
// Compiler:      microsoft visual c++ 6.0, unix/linux gcc
// System/OS:     Windows 32, Sun solaris, Linux, other unix systems (untested)
// Database:      none
// Description:   class CThermAlign - Thermodynamic Alignment Algorithm
// Author:        kaderali
// Date:          9/2000 - 12/2000
// Copyright:     (c) L. Kaderali, 9/2000 - 12/2000
//
// Revision History
// $ 00sep04 LK : created
// $ 00dec30 LK : changed to do local alignment of probe against
//                one entire sequence
// $ 01feb07 LK : optimized
// #$
//=============================================================================

#if !defined(AFX_THERMALIGN_H__1B9227F7_82AB_11D4_9FFF_000000000000__INCLUDED_)
#define AFX_THERMALIGN_H__1B9227F7_82AB_11D4_9FFF_000000000000__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <iostream>
#include <strstream>
//#include <ostream>
#include "nnparams.h"  // Nearest Neighbor Parameters
//#include "../gsfxtree/gsfxtree.h"             // Suffix Tree Stuff

using namespace std;

#ifdef _pack
#pragma pack(1)
#endif

//-----------------------------------------------------------------------------
// class CThermAlign

typedef class CThermAlign* PThermAlign;

class CThermAlign  
{
public:
	CThermAlign(int, int, PNNParams);
	virtual ~CThermAlign();
	void InitStrings(char*,char*,int,int,int=0);
	void InitBorder();
	//void CalculateCell(int, int);
	void CalculateTable(int=1);
	#ifdef _output_alignment
	    void printAlignment(ostream &outputStream);
		void OutputAlignment(ostream &outputStream);
		bool OutputAlignment(ostream &outputStream, int, int, int, bool local=false);
		void OutputLocalAlignment(ostream &outputStream);
	    void PrintDPTable(ostream&);
	#endif
	float GetEntropy(int, int);
	float GetEnthalpy(int, int);
	float GetFreeEnergy(int, int);
	float GetFreeEnergyK(int, int, float);
	float GetFreeEnergyC(int, int, float);
	float GetMeltingTempC(int, int);
	float GetMeltingTempK(int, int);
// lk01feb07: removed maxlocstuff as not requiered by thermtreealign...
	float maxloctm;      // maximum local temperature found
	int maxloci;         // i position thereof
	int maxlocj;         // j position thereof
	int maxloct;         // and type of maximum alignment!
    int targetNumber;    // id number of target sequence
//    PGSfxLeaf probeNode; // identifier for probe sequence
	char *seq1;         // Sequence 1
	char *seq2;         // Sequence 2
	int seq1len;        // Length of Sequence 1
	int seq2len;        // Length of Sequence 2
	int maxseq1len;     // length of longest target...
private:
	float *dH;         // Dynamic Programming Table for Entropy
	float *dS;         // Dynamic Programming Table for Enthalpy
	PNNParams NNParams; // Nearest Neighbor parameters
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

#endif // !defined(AFX_THERMALIGN_H__1B9227F7_82AB_11D4_9FFF_000000000000__INCLUDED_)
