//=============================================================================
// Module:        nnparams.cpp
// Project:       Diploma Thesis - Probe Selection for DNA Microarrays
// Type:          implementation - Nearest Neighbor Model / Parameters
// Language:      c++
// Compiler:      microsoft visual c++ 6.0, unix/linux gcc
// System/OS:     Windows 32, Sun solaris, Linux, other unix systems (untested)
// Database:      none
// Description:   class CNNParams - Nearest Neighbor Model / Parameters
// Author:        kaderali
// Date:          12/2000
// Copyright:     (c) L. Kaderali, 9/2000 - 12/2000
//
// Revision History
// $              00sep07 LK : created
//                00dec17 LK : modified to use new nn parameters
//                00dec29 LK : modified to include dangling end parameters
//                01jan09 LK : included CalcSelfTM
//                01feb07 LK : optimized
// #$
//=============================================================================

#include <memory.h>
#include <math.h>
#include <string.h>
#include"nnparams.h"

#ifdef _pack
#pragma pack(1)
#endif



float forbidden_entropy;

///////////////////////////////////////////////////////////////////////////////
// Construction/Destruction
///////////////////////////////////////////////////////////////////////////////

CNNParams::CNNParams()
{
	// initialize parameter table!
	InitParams();
}

CNNParams::~CNNParams()
{
	// does nothing!
}

///////////////////////////////////////////////////////////////////////////////
//  FUNCTION:   void CNNParams::InitParams()
//
//  PURPOSE:    Initialize nearest neighbor parameters. For now, simply set
//              params as requiered. Extend to read from file sometime.
//
//  PARAMETERS:
//      none
//
//  RETURN VALUE:
//      void
//
//  REVISION HISTORY
//  $              00sep06 : created LK
//                 00dec17 : modified to use new parameters LK
//                 00dec29 : included dangling end data LK
//  #$

void CNNParams::InitParams(float c1, float c2, float kp, int sm)
{
        Ct1 = c1;
	Ct2 = c2; 
	kplus = kp;
	int maxCT = 1;
	if(Ct2 > Ct1)
	{
		maxCT = 2;
	}
	float ctFactor;
 	if(Ct1 == Ct2) 
	{
		ctFactor = Ct1/2;
	}
	else if (maxCT == 1) 
	{
		ctFactor = Ct1-Ct2/2;
	}
	else
	{
		ctFactor = Ct2-Ct1/2;
	}
	rlogc = R * log(ctFactor);
	forbidden_entropy = rlogc;
	kfac = 0.368 * log (kplus);
	saltMethod = sm;
	int x,y,a,b;  // variables used as counters...

	// Set all parameters to zero!
	memset(dH,0,sizeof(dH));
	memset(dS,0,sizeof(dS));

	// Set all X-/Y-, -X/Y- and X-/-Y so, that TM will be VERY small! 
	for (x=1;x<=4;x++)
	{
		for (y=1;y<=4;y++)
		{
			ndH(0,x,y,0)=forbidden_enthalpy;
			ndS(0,x,y,0)=forbidden_entropy;
			ndH(x,0,0,y)=forbidden_enthalpy;
			ndS(x,0,0,y)=forbidden_entropy;
			ndH(x,0,y,0)=forbidden_enthalpy;
			ndS(x,0,y,0)=forbidden_entropy;
			// forbid X-/Y$ and X$/Y- etc., i.e. terminal must not be paired with gap!
			ndH(x,5,y,0)=forbidden_enthalpy;
			ndS(x,5,y,0)=forbidden_entropy;
			ndH(x,0,y,5)=forbidden_enthalpy;
			ndS(x,0,y,5)=forbidden_entropy;
			ndH(5,x,0,y)=forbidden_enthalpy;
			ndS(5,x,0,y)=forbidden_entropy;
			ndH(0,x,5,y)=forbidden_enthalpy;
			ndS(0,x,5,y)=forbidden_entropy;
			// forbid X$/-Y etc.
			ndH(x,5,0,y)=forbidden_enthalpy;
			ndS(x,5,0,y)=forbidden_entropy;
			ndH(x,0,5,y)=forbidden_enthalpy;
			ndS(x,0,5,y)=forbidden_entropy;
			ndH(5,x,y,0)=forbidden_enthalpy;
			ndS(5,x,y,0)=forbidden_entropy;
			ndH(0,x,y,5)=forbidden_enthalpy;
			ndS(0,x,y,5)=forbidden_entropy;
		
		}
		// also, forbid x-/-- and --/x-, i.e. no two inner gaps paired
		ndH(x,0,0,0)=forbidden_enthalpy;
		ndS(x,0,0,0)=forbidden_entropy;
		ndH(0,0,x,0)=forbidden_enthalpy;
		ndS(0,0,x,0)=forbidden_entropy;
		// x-/-$
		ndH(x,0,0,5)=forbidden_enthalpy;
		ndS(x,0,0,5)=forbidden_entropy;
		ndH(5,0,0,x)=forbidden_enthalpy;
		ndS(5,0,0,x)=forbidden_entropy;
		ndH(0,5,x,0)=forbidden_enthalpy;
		ndS(x,0,0,5)=forbidden_entropy;
		ndH(0,x,5,0)=forbidden_enthalpy;
		ndS(0,x,5,0)=forbidden_entropy;
	}
	// forbid --/--
	ndH(0,0,0,0)=forbidden_enthalpy;
	ndS(0,0,0,0)=forbidden_entropy;

	ndH(5,0,0,0)=forbidden_enthalpy;
	ndS(5,0,0,0)=forbidden_entropy;
	ndH(0,0,5,0)=forbidden_enthalpy;
	ndS(0,0,5,0)=forbidden_entropy;
	ndH(0,5,5,0)=forbidden_enthalpy;
	ndS(0,5,5,0)=forbidden_entropy;

	// Interior loops (double Mismatches)
    #define iloop_entropy -0.97f
    #define iloop_enthalpy 0.0f
	for (x=1; x<=4; x++)
		for (y=1; y<=4; y++)
			for (a=1; a<=4; a++)
				for (b=1; b<=4; b++)
					// AT and CG pair, and as A=1, C=2, G=3, T=4 this means
					// we have Watson-Crick pairs if (x+a==5) and (y+b)==5.
					if (!((x+a==5)||(y+b==5)))
					{
						// No watson-crick-pair, i.e. double mismatch!
						// set enthalpy/entropy to loop expansion!
						ndH(x,y,a,b) = iloop_enthalpy;
						ndS(x,y,a,b) = iloop_entropy;
					}

	// xy/-- and --/xy (Bulge Loops of size > 1)
    #define bloop_entropy -1.3f
    #define bloop_enthalpy 0.0f
	for (x=1; x<=4; x++)
		for (y=1; y<=4; y++)
		{
			ndH(x,y,0,0) = bloop_enthalpy;
			ndS(x,y,0,0) = bloop_entropy;
			ndH(0,0,x,y) = bloop_enthalpy;
			ndS(0,0,x,y) = bloop_entropy;
		}

    // x-/ya abd xa/y- as well as -x/ay and ax/-y
	// bulge opening and closing parameters with
	// adjacent matches / mismatches
	// obulge_mism and cbulge_mism chosen so high to avoid 
	//     AAAAAAAAA
	//     T--G----T
	// being better than
	//     AAAAAAAAA
	//     TG------T
    #define obulge_match_H (-2.66f * 1000)
    #define obulge_match_S -14.22f
    #define cbulge_match_H (-2.66f * 1000)
    #define cbulge_match_S -14.22f
    #define obulge_mism_H (0.0f * 1000)
    #define obulge_mism_S -6.45f
    #define cbulge_mism_H 0.0f
    #define cbulge_mism_S -6.45f
	for (x=1; x<=4; x++)
		for (y=1; y<=4; y++)
			for (a=1; a<=4; a++)
			{
				if (x+y==5)  // other base pair matches!
				{
					ndH(x,0,y,a)=obulge_match_H;  // bulge opening
					ndS(x,0,y,a)=obulge_match_S;
					ndH(x,a,y,0)=obulge_match_H;
					ndS(x,a,y,0)=obulge_match_S;
					ndH(0,x,a,y)=cbulge_match_H;  // bulge closing
					ndS(0,x,a,y)=cbulge_match_S;
					ndH(a,x,0,y)=cbulge_match_H;
					ndS(a,x,0,y)=cbulge_match_S;
				}
				else
				{           // mismatch in other base pair!
					ndH(x,0,y,a)=obulge_mism_H;   // bulge opening
					ndS(x,0,y,a)=obulge_mism_S;
					ndH(x,a,y,0)=obulge_mism_H;
					ndS(x,a,y,0)=obulge_mism_S;
					ndH(0,x,a,y)=cbulge_mism_H;   // bulge closing
					ndS(0,x,a,y)=cbulge_mism_S;
					ndH(a,x,0,y)=cbulge_mism_H;
					ndS(a,x,0,y)=cbulge_mism_S;
				}
			}

	// Watson-Crick pairs (note that only ten are unique, as obviously
	// 5'-AG-3'/3'-TC-5'  =  5'-CT-3'/3'-GA-5' etc.
	ndH(1,1,4,4)=-7.6f*1000;  ndS(1,1,4,4)=-21.3f;   // AA/TT 04
	ndH(1,2,4,3)=-8.4f*1000;  ndS(1,2,4,3)=-22.4f;   // AC/TG adapted GT/CA
	ndH(1,3,4,2)=-7.8f*1000;  ndS(1,3,4,2)=-21.0f;   // AG/TC adapted CT/GA
	ndH(1,4,4,1)=-7.2f*1000;  ndS(1,4,4,1)=-20.4f;   // AT/TA 04
	ndH(2,1,3,4)=-8.5f*1000;  ndS(2,1,3,4)=-22.7f;   // CA/GT 04
	ndH(2,2,3,3)=-8.0f*1000;  ndS(2,2,3,3)=-19.9f;   // CC/GG adapted GG/CC
	ndH(2,3,3,2)=-10.6f*1000; ndS(2,3,3,2)=-27.2f;   // CG/GC 04
	ndH(2,4,3,1)=-7.8f*1000;  ndS(2,4,3,1)=-21.0f;   // CT/GA 04
	ndH(3,1,2,4)=-8.2f*1000;  ndS(3,1,2,4)=-22.2f;   // GA/CT 04
	ndH(3,2,2,3)=-9.8f*1000;  ndS(3,2,2,3)=-24.4f;   // GC/CG 04
	ndH(3,3,2,2)=-8.0f*1000;  ndS(3,3,2,2)=-19.9f;   // GG/CC 04
	ndH(3,4,2,1)=-8.4f*1000;  ndS(3,4,2,1)=-22.4f;   // GT/CA 04
	ndH(4,1,1,4)=-7.2f*1000;  ndS(4,1,1,4)=-21.3f;   // TA/AT 04
	ndH(4,2,1,3)=-8.2f*1000;  ndS(4,2,1,3)=-22.2f;   // TC/AG adapted GA/CT
	ndH(4,3,1,2)=-8.5f*1000;  ndS(4,3,1,2)=-22.7f;   // TG/AC adapted CA/GT
	ndH(4,4,1,1)=-7.6f*1000;  ndS(4,4,1,1)=-21.3f;   // TT/AA adapted AA/TT

    // A-C Mismatches (Values for pH 7.0)
	ndH(1,1,2,4)=7.6f*1000;   ndS(1,1,2,4)=20.2f;    // AA/CT
	ndH(1,1,4,2)=2.3f*1000;   ndS(1,1,4,2)=4.6f;     // AA/TC
	ndH(1,2,2,3)=-0.7f*1000;  ndS(1,2,2,3)=-3.8f;    // AC/CG
	ndH(1,2,4,1)=5.3f*1000;   ndS(1,2,4,1)=14.6f;    // AC/TA
	ndH(1,3,2,2)=0.6f*1000;   ndS(1,3,2,2)=-0.6f;    // AG/CC
	ndH(1,4,2,1)=5.3f*1000;   ndS(1,4,2,1)=14.6f;    // AT/CA
	ndH(2,1,1,4)=3.4f*1000;   ndS(2,1,1,4)=8.0f;     // CA/AT
	ndH(2,1,3,2)=1.9f*1000;   ndS(2,1,3,2)=3.7f;     // CA/GC
	ndH(2,2,1,3)=5.2f*1000;   ndS(2,2,1,3)=14.2f;    // CC/AG
	ndH(2,2,3,1)=0.6f*1000;   ndS(2,2,3,1)=-0.6f;    // CC/GA
	ndH(2,3,1,2)=1.9f*1000;   ndS(2,3,1,2)=3.7f;     // CG/AC
	ndH(2,4,1,1)=2.3f*1000;   ndS(2,4,1,1)=4.6f;     // CT/AA
	ndH(3,1,2,2)=5.2f*1000;   ndS(3,1,2,2)=14.2f;    // GA/CC
	ndH(3,2,2,1)=-0.7f*1000;  ndS(3,2,2,1)=-3.8f;    // GC/CA
	ndH(4,1,1,2)=3.4f*1000;   ndS(4,1,1,2)=8.0f;     // TA/AC
	ndH(4,2,1,1)=7.6f*1000;   ndS(4,2,1,1)=20.2f;    // TC/AA

	// C-T Mismatches
	ndH(1,2,4,4)=0.7f*1000;   ndS(1,2,4,4)=0.2f;     // AC/TT
	ndH(1,4,4,2)=-1.2f*1000;  ndS(1,4,4,2)=-6.2f;    // AT/TC
	ndH(2,1,4,4)=1.0f*1000;   ndS(2,1,4,4)=0.7f;     // CA/TT
	ndH(2,2,3,4)=-0.8f*1000;  ndS(2,2,3,4)=-4.5f;    // CC/GT
	ndH(2,2,4,3)=5.2f*1000;   ndS(2,2,4,3)=13.5f;    // CC/TG
	ndH(2,3,4,2)=-1.5f*1000;  ndS(2,3,4,2)=-6.1f;    // CG/TC
	ndH(2,4,3,2)=-1.5f*1000;  ndS(2,4,3,2)=-6.1f;    // CT/GC
	ndH(2,4,4,1)=-1.2f*1000;  ndS(2,4,4,1)=-6.2f;    // CT/TA
	ndH(3,2,2,4)=2.3f*1000;   ndS(3,2,2,4)=5.4f;     // GC/CT
	ndH(3,4,2,2)=5.2f*1000;   ndS(3,4,2,2)=13.5f;    // GT/CC
	ndH(4,1,2,4)=1.2f*1000;   ndS(4,1,2,4)=0.7f;     // TA/CT
	ndH(4,2,2,3)=2.3f*1000;   ndS(4,2,2,3)=5.4f;     // TC/CG
	ndH(4,2,1,4)=1.2f*1000;   ndS(4,2,1,4)=0.7f;     // TC/AT
	ndH(4,3,2,2)=-0.8f*1000;  ndS(4,3,2,2)=-4.5f;    // TG/CC
	ndH(4,4,2,1)=0.7f*1000;   ndS(4,4,2,1)=0.2f;     // TT/CA
	ndH(4,4,1,2)=1.0f*1000;   ndS(4,4,1,2)=0.7f;     // TT/AC

	// G-A Mismatches
	ndH(1,1,3,4)=3.0f*1000;   ndS(1,1,3,4)=7.4f;     // AA/GT
	ndH(1,1,4,3)=-0.6f*1000;  ndS(1,1,4,3)=-2.3f;    // AA/TG
	ndH(1,2,3,3)=0.5f*1000;   ndS(1,2,3,3)=3.2f;     // AC/GG
	ndH(1,3,3,2)=-4.0f*1000;  ndS(1,3,3,2)=-13.2f;   // AG/GC
	ndH(1,3,4,1)=-0.7f*1000;  ndS(1,3,4,1)=-2.3f;    // AG/TA
	ndH(1,4,3,1)=-0.7f*1000;  ndS(1,4,3,1)=-2.3f;    // AT/GA
	ndH(2,1,3,3)=-0.7f*1000;  ndS(2,1,3,3)=-2.3f;    // CA/GG
	ndH(2,3,3,1)=-4.0f*1000;  ndS(2,3,3,1)=-13.2f;   // CG/GA
	ndH(3,1,1,4)=0.7f*1000;   ndS(3,1,1,4)=0.7f;     // GA/AT
	ndH(3,1,2,3)=-0.6f*1000;  ndS(3,1,2,3)=-1.0f;    // GA/CG
	ndH(3,2,1,3)=-0.6f*1000;  ndS(3,2,1,3)=-1.0f;    // GC/AG
	ndH(3,3,1,2)=-0.7f*1000;  ndS(3,3,1,2)=-2.3f;    // GG/AC
	ndH(3,3,2,1)=0.5f*1000;   ndS(3,3,2,1)=3.2f;     // GG/CA
	ndH(3,4,1,1)=-0.6f*1000;  ndS(3,4,1,1)=-2.3f;    // GT/AA
	ndH(4,1,1,3)=0.7f*1000;   ndS(4,1,1,3)=0.7f;     // TA/AG
	ndH(4,3,1,1)=3.0f*1000;   ndS(4,3,1,1)=7.4f;     // TG/AA

	// G-T Mismatches
	ndH(1,3,4,4)=1.0f*1000;   ndS(1,3,4,4)=0.9f;     // AG/TT
	ndH(1,4,4,3)=-2.5f*1000;  ndS(1,4,4,3)=-8.3f;    // AT/TG
	ndH(2,3,3,4)=-4.1f*1000;  ndS(2,3,3,4)=-11.7f;   // CG/GT
	ndH(2,4,3,3)=-2.8f*1000;  ndS(2,4,3,3)=-8.0f;    // CT/GG
	ndH(3,1,4,4)=-1.3f*1000;  ndS(3,1,4,4)=-5.3f;    // GA/TT
	ndH(3,2,4,3)=-4.4f*1000;  ndS(3,2,4,3)=-12.3f;   // GC/TG
	ndH(3,3,2,4)=3.3f*1000;   ndS(3,3,2,4)=10.4f;    // GG/CT
	ndH(3,3,4,2)=-2.8f*1000;  ndS(3,3,4,2)=-8.0f;    // GG/TC
//	ndH(3,3,4,4)=5.8f*1000;   ndS(3,3,4,4)=16.3f;    // GG/TT
	ndH(3,4,2,3)=-4.4f*1000;  ndS(3,4,2,3)=-12.3f;   // GT/CG
	ndH(3,4,4,1)=-2.5f*1000;  ndS(3,4,4,1)=-8.3f;    // GT/TA
//	ndH(3,4,4,3)=4.1f*1000;   ndS(3,4,4,3)=9.5f;     // GT/TG
	ndH(4,1,3,4)=-0.1f*1000;  ndS(4,1,3,4)=-1.7f;    // TA/GT
	ndH(4,2,3,3)=3.3f*1000;   ndS(4,2,3,3)=10.4f;    // TC/GG
	ndH(4,3,1,4)=-0.1f*1000;  ndS(4,3,1,4)=-1.7f;    // TG/AT
	ndH(4,3,3,2)=-4.1f*1000;  ndS(4,3,3,2)=-11.7f;   // TG/GC
//	ndH(4,3,3,4)=-1.4f*1000;  ndS(4,3,3,4)=-6.2f;    // TG/GT
	ndH(4,4,1,3)=-1.3f*1000;  ndS(4,4,1,3)=-5.3f;    // TT/AG
	ndH(4,4,3,1)=1.0f*1000;   ndS(4,4,3,1)=0.9f;     // TT/GA
//	ndH(4,4,3,3)=5.8f*1000;   ndS(4,4,3,3)=16.3f;    // TT/GG

	// A-A Mismatches
	ndH(1,1,1,4)=4.7f*1000;   ndS(1,1,1,4)=12.9f;    // AA/AT
	ndH(1,1,4,1)=1.2f*1000;   ndS(1,1,4,1)=1.7f;     // AA/TA
	ndH(1,2,1,3)=-2.9f*1000;  ndS(1,2,1,3)=-9.8f;    // AC/AG
	ndH(1,3,1,2)=-0.9f*1000;  ndS(1,3,1,2)=-4.2f;    // AG/AC
	ndH(1,4,1,1)=1.2f*1000;   ndS(1,4,1,1)=1.7f;     // AT/AA
	ndH(2,1,3,1)=-0.9f*1000;  ndS(2,1,3,1)=-4.2f;    // CA/GA
    ndH(3,1,2,1)=-2.9f*1000;  ndS(3,1,2,1)=-9.8f;    // GA/CA
	ndH(4,1,1,1)=4.7f*1000;   ndS(4,1,1,1)=12.9f;    // TA/AA

	// C-C Mismatches
	ndH(1,2,4,2)=0.0f*1000;   ndS(1,2,4,2)=-4.4f;    // AC/TC
	ndH(2,1,2,4)=6.1f*1000;   ndS(2,1,2,4)=16.4f;    // CA/CT
	ndH(2,2,2,3)=3.6f*1000;   ndS(2,2,2,3)=8.9f;     // CC/CG
	ndH(2,2,3,2)=-1.5f*1000;  ndS(2,2,3,2)=-7.2f;    // CC/GC
	ndH(2,3,2,2)=-1.5f*1000;  ndS(2,3,2,2)=-7.2f;    // CG/CC
	ndH(2,4,2,1)=0.0f*1000;   ndS(2,4,2,1)=-4.4f;    // CT/CA
	ndH(3,2,2,2)=3.6f*1000;   ndS(3,2,2,2)=8.9f;     // GC/CC
	ndH(4,2,1,2)=6.1f*1000;   ndS(4,2,1,2)=16.4f;    // TC/AC

	// G-G Mismatches
	ndH(1,3,4,3)=-3.1f*1000;  ndS(1,3,4,3)=-9.5f;    // AG/TG
	ndH(2,3,3,3)=-4.9f*1000;  ndS(2,3,3,3)=-15.3f;   // CG/GG
	ndH(3,1,3,4)=1.6f*1000;   ndS(3,1,3,4)=3.6f;     // GA/GT
	ndH(3,2,3,3)=-6.0f*1000;  ndS(3,2,3,3)=-15.8f;   // GC/GG
	ndH(3,3,2,3)=-6.0f*1000;  ndS(3,3,2,3)=-15.8f;   // GG/CG
	ndH(3,3,3,2)=-4.9f*1000;  ndS(3,3,3,2)=-15.3f;   // GG/GC
	ndH(3,4,3,1)=-3.1f*1000;  ndS(3,4,3,1)=-9.5f;    // GT/GA
	ndH(4,3,1,3)=1.6f*1000;   ndS(4,3,1,3)=3.6f;     // TG/AG

	// T-T Mismatches
	ndH(1,4,4,4)=-2.7f*1000;  ndS(1,4,4,4)=-10.8f;   // AT/TT
	ndH(2,4,3,4)=-5.0f*1000;  ndS(2,4,3,4)=-15.8f;   // CT/GT
	ndH(3,4,2,4)=-2.2f*1000;  ndS(3,4,2,4)=-8.4f;    // GT/CT
	ndH(4,1,4,4)=0.2f*1000;   ndS(4,1,4,4)=-1.5f;    // TA/TT
	ndH(4,2,4,3)=-2.2f*1000;  ndS(4,2,4,3)=-8.4f;    // TC/TG
	ndH(4,3,4,2)=-5.0f*1000;  ndS(4,3,4,2)=-15.8f;   // TG/TC
	ndH(4,4,1,4)=0.2f*1000;   ndS(4,4,1,4)=-1.5f;    // TT/AT
	ndH(4,4,4,1)=-2.7f*1000;  ndS(4,4,4,1)=-10.8f;   // TT/TA

	// Dangling Ends
	ndH(5,1,1,4)=-0.7f*1000;  ndS(5,1,1,4)=-0.8f;    // $A/AT
	ndH(5,1,2,4)=4.4f*1000;   ndS(5,1,2,4)=14.9f;    // $A/CT
	ndH(5,1,3,4)=-1.6f*1000;  ndS(5,1,3,4)=-3.6f;    // $A/GT
	ndH(5,1,4,4)=2.9f*1000;   ndS(5,1,4,4)=10.4f;    // $A/TT
	ndH(5,2,1,3)=-2.1f*1000;  ndS(5,2,1,3)=-3.9f;    // $C/AG
	ndH(5,2,2,3)=-0.2f*1000;  ndS(5,2,2,3)=-0.1f;    // $C/CG
	ndH(5,2,3,3)=-3.9f*1000;  ndS(5,2,3,3)=-11.2f;   // $C/GG
	ndH(5,2,4,3)=-4.4f*1000;  ndS(5,2,4,3)=-13.1f;   // $C/TG
	ndH(5,3,1,2)=-5.9f*1000;  ndS(5,3,1,2)=-16.5f;   // $G/AC
	ndH(5,3,2,2)=-2.6f*1000;  ndS(5,3,2,2)=-7.4f;    // $G/CC
	ndH(5,3,3,2)=-3.2f*1000;  ndS(5,3,3,2)=-10.4f;   // $G/GC
	ndH(5,3,4,2)=-5.2f*1000;  ndS(5,3,4,2)=-15.0f;   // $G/TC
	ndH(5,4,1,1)=-0.5f*1000;  ndS(5,4,1,1)=-1.1f;    // $T/AA
	ndH(5,4,2,1)=4.7f*1000;   ndS(5,4,2,1)=14.2f;    // $T/CA
	ndH(5,4,3,1)=-4.1f*1000;  ndS(5,4,3,1)=-13.1f;   // $T/GA
	ndH(5,4,4,1)=-3.8f*1000;  ndS(5,4,4,1)=-12.6f;   // $T/TA
	ndH(1,5,4,1)=-2.9f*1000;  ndS(1,5,4,1)=-7.6f;    // A$/TA
	ndH(1,5,4,2)=-4.1f*1000;  ndS(1,5,4,2)=-13.0f;   // A$/TC
	ndH(1,5,4,3)=-4.2f*1000;  ndS(1,5,4,3)=-15.0f;   // A$/TG
	ndH(1,5,4,4)=-0.2f*1000;  ndS(1,5,4,4)=-0.5f;    // A$/TT
	ndH(1,1,5,4)=0.2f*1000;   ndS(1,1,5,4)=2.3f;     // AA/$T
	ndH(1,1,4,5)=-0.5f*1000;  ndS(1,1,4,5)=-1.1f;    // AA/T$
	ndH(1,2,5,3)=-6.3f*1000;  ndS(1,2,5,3)=-17.1f;   // AC/$G
	ndH(1,2,4,5)=4.7f*1000;   ndS(1,2,4,5)=14.2f;    // AC/T$
	ndH(1,3,5,2)=-3.7f*1000;  ndS(1,3,5,2)=-10.0f;   // AG/$C
	ndH(1,3,4,5)=-4.1f*1000;  ndS(1,3,4,5)=-13.1f;   // AG/T$
	ndH(1,4,5,1)=-2.9f*1000;  ndS(1,4,5,1)=-7.6f;    // AT/$A
	ndH(1,4,4,5)=-3.8f*1000;  ndS(1,4,4,5)=-12.6f;   // AT/T$
	ndH(2,5,3,1)=-3.7f*1000;  ndS(2,5,3,1)=-10.0f;   // C$/GA
	ndH(2,5,3,2)=-4.0f*1000;  ndS(2,5,3,2)=-11.9f;   // C$/GC
	ndH(2,5,3,3)=-3.9f*1000;  ndS(2,5,3,3)=-10.9f;   // C$/GG
	ndH(2,5,3,4)=-4.9f*1000;  ndS(2,5,3,4)=-13.8f;   // C$/GT
	ndH(2,1,5,4)=0.6f*1000;   ndS(2,1,5,4)=3.3f;     // CA/$T
	ndH(2,1,3,5)=-5.9f*1000;  ndS(2,1,3,5)=-16.5f;   // CA/G$
	ndH(2,2,5,3)=-4.4f*1000;  ndS(2,2,5,3)=-12.6f;   // CC/$G
	ndH(2,2,3,5)=-2.6f*1000;  ndS(2,2,3,5)=-7.4f;    // CC/G$
	ndH(2,3,5,2)=-4.0f*1000;  ndS(2,3,5,2)=-11.9f;   // CG/$C
	ndH(2,3,3,5)=-3.2f*1000;  ndS(2,3,3,5)=-10.4f;   // CG/G$
	ndH(2,4,5,1)=-4.1f*1000;  ndS(2,4,5,1)=-13.0f;   // CT/$A
	ndH(2,4,3,5)=-5.2f*1000;  ndS(2,4,3,5)=-15.0f;   // CT/G$
	ndH(3,5,2,1)=-6.3f*1000;  ndS(3,5,2,1)=-17.1f;   // G$/CA
	ndH(3,5,2,2)=-4.4f*1000;  ndS(3,5,2,2)=-12.6f;   // G$/CC
	ndH(3,5,2,3)=-5.1f*1000;  ndS(3,5,2,3)=-14.0f;   // G$/CG
	ndH(3,5,2,4)=-4.0f*1000;  ndS(3,5,2,4)=-10.9f;   // G$/CT
	ndH(3,1,5,4)=-1.1f*1000;  ndS(3,1,5,4)=-1.6f;    // GA/$T
	ndH(3,1,2,5)=-2.1f*1000;  ndS(3,1,2,5)=-3.9f;    // GA/C$
	ndH(3,2,5,3)=-5.1f*1000;  ndS(3,2,5,3)=-14.0f;   // GC/$G
	ndH(3,2,2,5)=-0.2f*1000;  ndS(3,2,2,5)=-0.1f;    // GC/C$
	ndH(3,3,5,2)=-3.9f*1000;  ndS(3,3,5,2)=-10.9f;   // GG/$C
	ndH(3,3,2,5)=-3.9f*1000;  ndS(3,3,2,5)=-11.2f;   // GG/C$
	ndH(3,4,5,1)=-4.2f*1000;  ndS(3,4,5,1)=-15.0f;   // GT/$A
	ndH(3,4,2,5)=-4.4f*1000;  ndS(3,4,2,5)=-13.1f;   // GT/C$
	ndH(4,5,1,1)=0.2f*1000;   ndS(4,5,1,1)=2.3f;     // T$/AA
	ndH(4,5,1,2)=0.6f*1000;   ndS(4,5,1,2)=3.3f;     // T$/AC
	ndH(4,5,1,3)=-1.1f*1000;  ndS(4,5,1,3)=-1.6f;    // T$/AG
	ndH(4,5,1,4)=-6.9f*1000;  ndS(4,5,1,4)=-20.0f;   // T$/AT
	ndH(4,1,5,4)=-6.9f*1000;  ndS(4,1,5,4)=-20.0f;   // TA/$T
	ndH(4,1,1,5)=-0.7f*1000;  ndS(4,1,1,5)=-0.7f;    // TA/A$
	ndH(4,2,5,3)=-4.0f*1000;  ndS(4,2,5,3)=-10.9f;   // TC/$G
	ndH(4,2,1,5)=4.4f*1000;   ndS(4,2,1,5)=14.9f;    // TC/A$
	ndH(4,3,5,2)=-4.9f*1000;  ndS(4,3,5,2)=-13.8f;   // TG/$C
	ndH(4,3,1,5)=-1.6f*1000;  ndS(4,3,1,5)=-3.6f;    // TG/A$
	ndH(4,4,5,1)=-0.2f*1000;  ndS(4,4,5,1)=-0.5f;    // TT/$A
	ndH(4,4,1,5)=2.9f*1000;   ndS(4,4,1,5)=10.4f;    // TT/A$

	return;
}

void CNNParams::UpdateParams(char * s1, char * s2)
{
	float l1 = strlen(s1);
	float l2 = strlen(s2);
	if(l1<l2) {
		gcContent = CountGCContent(s1);
		gcContent /= (l1-2);
	} else if (l1>l2) {
		gcContent = CountGCContent(s2);
		gcContent /= (l2-2);
	} else {
		gcContent = CountGCContent(s1)+CountGCContent(s2);
		gcContent /= (l1+l2-4);
	}
}

float CNNParams::CountGCContent(char * seq ) {
	int lseq = strlen(seq);
	float count = 0;
	for(int k=0;k<lseq;k++) {
		if (seq[k] == 'G' || seq[k] == 'C' ) {
			count+=1;
		}
	}
	return count;
	
}

///////////////////////////////////////////////////////////////////////////////
//  FUNCTION:   static char CNNParams::getComplement(char mychar)
//
//  PURPOSE:    return Watson-Crick complement to mychar
//
//  PARAMETERS:
//	    mychar - Character to complement
//
//  RETURN VALUE:
//      char   - complement to mychar
//
//  REVISION HISTORY
//  $              00sep28 : created LK
//  #$

char CNNParams::getComplement(char mychar, bool asnum)
{
	if (mychar==1)  // A -> T
		return 4;
	if (mychar==2)  // C -> G
		return 3;
	if (mychar==3)  // G -> T
		return 2;
	if (mychar==4)  // T -> A
		return 1;
	if (mychar==5)  // $ -> $
		return 5;
	if (!asnum)
	{
		if (mychar=='A')  // A -> T
			return 'T';
		if (mychar=='C')  // C -> G
			return 'G';
		if (mychar=='G')  // G -> T
			return 'C';
		if (mychar=='T')  // T -> A
			return 'A';
		if (mychar=='$')  // $ -> $
			return '$';
	}
	else
	{
		if (mychar=='A')  // A -> T
			return 4;
		if (mychar=='C')  // C -> G
			return 3;
		if (mychar=='G')  // G -> T
			return 2;
		if (mychar=='T')  // T -> A
			return 1;
		if (mychar=='$')  // $ -> $
			return 5;
	}
	return '*';
}



///////////////////////////////////////////////////////////////////////////////
//  FUNCTION:   static bool CNNParams::isMismatch(char,char)
//
//  PURPOSE:    Return true if char1 and char2 are not watson-crick pair
//
//  PARAMETERS:
//      char1   - first character
//      char2   - second character
//
//  RETURN VALUE:
//      bool    - true if char1,char2 are Watson-Crick, false otw
//
//  REVISION HISTORY
//  $              00sep28 : created LK
//  #$

bool CNNParams::isMismatch(char c1,char c2)
{
	if (c1+c2==5)
		return false;
	else
		return true;
}




