//=============================================================================
// Module:        nnparams.h
// Project:       Diploma Thesis - Probe Selection for DNA Microarrays
// Type:          header file - Nearest Neighbor Parameters / Model.
// Language:      c++
// Compiler:      microsoft visual c++ 6.0, unix/linux gcc
// System/OS:     Windows 32, Sun solaris, Linux, other unix systems (untested)
// Database:      none
// Description:   class CNNParams - Nearest Neighbor Model Parameters
// Author:        kaderali
// Date:          9-12/2000
// Copyright:     (c) L. Kaderali, 9/2000 - 12/2000
//
// Revision History
// $              00sep07 LK : created
//                00dec29 LK : changed to include dangling end data 
//                01jan09 LK : included CalcSelfTM function 
//                01feb07 LK : optimized
// #$
//=============================================================================

#if !defined(AFX_NNPARAMS_H__05604705_84E8_11D4_A001_000000000000__INCLUDED_)
#define AFX_NNPARAMS_H__05604705_84E8_11D4_A001_000000000000__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <math.h>

#ifdef _pack
#pragma pack(1)
#endif


// following defines to simplify coding...
#define ndH(a,b,c,d) dH[a][b][c][d]
#define ndS(a,b,c,d) dS[a][b][c][d]
#define forbidden_enthalpy 1000000000000000000.0f
//#define forbidden_enthalpy_div1000 1000000000000000.0f
// forbidden entropy=-rlogc
// #define forbidden_entropy 30.205986374220235304486574573422f
// Boltzmann factor (cal/degrees C*mol)
#define	R 1.987f
#define	SALT_METHOD_SANTALUCIA 1
#define	SALT_METHOD_OWCZARZY 2
// Strand concentration (assumption!) (M)
// #define Ct  0.000001f
// r*ln(ct/4) as required by many formulas
//#define rlogc -30.205986374220235304486574573422f
extern float forbidden_entropy;

//-----------------------------------------------------------------------------
// class CNNParams
typedef class CNNParams* PNNParams;


class CNNParams  
{
public:
        float Ct1;
	float Ct2;
	float rlogc;
	float kplus;
	float kfac;
	int saltMethod;
	float gcContent;
        float new_TM; // geändert von ML!!!
	CNNParams();
	virtual ~CNNParams();
	void InitParams(float c1=0.000001f, float c2=0.000001f, float kp=1, int sm = SALT_METHOD_SANTALUCIA );
	void UpdateParams(char * s1, char * s2);
	float CountGCContent(char * seq );

	inline char convertNum(char c)
	{
		if (c=='A')
			return 1;
		if (c=='C')
			return 2;
		if (c=='G')
			return 3;
		if (c=='T')
			return 4;
		if (c=='$')
          		return 5;
		return 0;
	}

	///////////////////////////////////////////////////////////////////////////////
	// inline function
	///////////////////////////////////////////////////////////////////////////////
	//  FUNCTION:   float CNNParams::GetEntropy(char,char,char,char)
	//
	//  PURPOSE:    Retrieve Entropy for given NN-Pair from parameter table
	//
	//  PARAMETERS: x0,x1,y0,y1: Pairs to look up in form  
	//              5'-x0-x1-3' / 3'-y0-y1-5'
	//
	//  RETURN VALUE: 
	//      float: Entropy dS for given NN pair
	//
	//  REVISION HISTORY
	//  $              00sep06 : created LK
	//                 00dec29 : included dangling end parameters
	//                 01feb07 : rewritten. looks ugly now, but is FAST! inline.
	//  #$

	inline float GetEntropy(char x0, char x1, char y0, char y1)
	{
		char nx0=convertNum(x0);
		char nx1=convertNum(x1);
		char ny0=convertNum(y0);
		char ny1=convertNum(y1);
		float answer = ndS(nx0,nx1,ny0,ny1);
		//Salt correction Santalucia
		if (saltMethod == SALT_METHOD_SANTALUCIA) {
			if(nx0!=5 && 1<= nx1 && nx1<=4) {
				answer += 0.5*kfac;
			}
			if(ny1!=5 && 1<= ny0 && ny0<=4) {
				answer += 0.5*kfac;
			}
		}
		//Salt correction Owczarzy
		if (saltMethod == SALT_METHOD_OWCZARZY) {
			float logk = log(kplus);
			answer += ndH(nx0,nx1,ny0,ny1)*((4.29 * gcContent-3.95)*0.00001*logk+ 0.0000094*logk*logk);
		}	
		return answer;
	}	
	
	
	///////////////////////////////////////////////////////////////////////////////
	// inline function
	///////////////////////////////////////////////////////////////////////////////
	//  FUNCTION:   float CNNParams::GetEnthalpy(char,char,char,char)
	//
	//  PURPOSE:    Retrieve Enthalpy for given NN-Pair from parameter table
	//
	//  PARAMETERS: x0,x1,y0,y1: Pairs to look up in form  
	//              5'-x0-x1-3' / 3'-y0-y1-5'
	//
	//  RETURN VALUE: 
	//      float: Enthalpy dH for given NN pair
	//
	//  REVISION HISTORY
	//  $ 00sep06 : created LK
	//  $ 00dec29 : included dangling end parameters
	//  $ 01feb07 : rewritten. looks ugly now, but is FAST! inline.
	//  #$

	inline float GetEnthalpy(char x0, char x1, char y0, char y1)
	{
		char nx0=convertNum(x0);
		char nx1=convertNum(x1);
		char ny0=convertNum(y0);
		char ny1=convertNum(y1);
		return ndH(nx0,nx1,ny0,ny1);
	}



    ///////////////////////////////////////////////////////////////////////////////
	// inline function
	///////////////////////////////////////////////////////////////////////////////
	//  FUNCTION:   float CNNParams::CalcTM(float entropy,float enthalpy)
	//
	//  PURPOSE:    Return melting temperature TM for given entropy and enthalpy
	//              Assuming a one-state transition and using the formula
	//                TM = dH / (dS + R ln(Ct/4))
	//              entropy = dS + R ln Ct/4 (must already be included!)
	//              enthaklpy = dH
	//              where
	//                dH = enthalpy
	//                dS = entropy
	//                R  = Boltzmann factor
	//                Ct = Strand Concentration
	//
	//  PARAMETERS:
	//
	//  RETURN VALUE:
	//
	//  REVISION HISTORY
	//  $ 00sep06 : created LK
	//  $ 01jan07 : modified and corrected
	//  $ 01feb07 : optimized!!! inline
	//  $ 01feb09 : changed to include r ln ct in entropy!!!
	//  #$

	inline float CalcTM(float entropy,float enthalpy)
	{
		float tm = 0;               // absolute zero - return if model fails!
	    if (enthalpy>=forbidden_enthalpy)   //||(entropy==-cfact))  
			return 0;
		if (entropy<0)         // avoid division by zero and model errors!
		{
		  tm = enthalpy/entropy;// - kfac; //LKFEB
			if (tm<0)
				return 0;
		}
		return tm;
	}



    ///////////////////////////////////////////////////////////////////////////////
	// inline function
	///////////////////////////////////////////////////////////////////////////////
	//  FUNCTION:   void CNNParams::AlterTM(float)
	//
	//  PURPOSE:    TM can be altered by a new assignment
        //
	//  PARAMETERS: tm_new: new value for TM
	//
	//  RETURN VALUE:
	//
	//  REVISION HISTORY
	//
	//  #$

        inline void AlterTM(float tm_new)
        {
            new_TM = tm_new;   
            return;
        }



    ///////////////////////////////////////////////////////////////////////////////
	// inline function
	///////////////////////////////////////////////////////////////////////////////
	//  FUNCTION:   float CNNParams::CalcG(float entropy, float enthalpy)
	//
	//  PURPOSE:    return free energy G for given entropy and enthalpy
        //              Assuming a one-state transition and using the formula
        //              G = dH - new_TM * dS
        //                dH = enthalpy
	//                dS = entropy
	//                new_TM = value for the optimal melting temperature of
  	//                         the last iteration
 	// 
	//  PARAMETERS: entropy and enthalpy
	//
	//  RETURN VALUE: free Energy value (dG) for the 
	//
	//  REVISION HISTORY
	//
	//  #$

        inline float CalcG(float entropy, float enthalpy)
        {
            float freeEnergy = -999999999;
	    if (enthalpy>=forbidden_enthalpy)  
			return -999999999;
                if (entropy<0)         // avoid division by zero and model errors! 
                {
                    entropy = entropy * -1;
                    enthalpy = enthalpy * -1;
                    freeEnergy = enthalpy - new_TM * entropy;
                }
                return freeEnergy;
        }




    ///////////////////////////////////////////////////////////////////////////////
	// inline function
	///////////////////////////////////////////////////////////////////////////////
	//  FUNCTION:   float CNNParams::CalcSelfTM(char*)
	//
	//  PURPOSE:    Calculate TM for given sequence against its complement
	//
	//  PARAMETERS:
	//      char*   - Sequence to consider
	//
	//  RETURN VALUE:
	//      float   - Melting temperature in degrees Kelvin
	//
	//  REVISION HISTORY
	//  $ 01jan09 LK : created
	//  $ 01feb07 LK : inline.
	//  #$

	inline float CalcSelfTM(char* seq)
	{
		float thedH = 0;
		float thedS = GetInitialEntropy();
		char c1;
		char c2;
		char c3;
		char c4;
		for (unsigned int i=1;i<strlen(seq);i++)
		{
			c1 = getComplement(seq[i-1],true);
			c2  = getComplement(seq[i],true);
			c3 = seq[i-1];
			c4 = seq[i];
			if (c3>5)
			{
				if (c3=='A')
					c3=1;
				if (c3=='C')
					c3=2;
				if (c3=='G')
					c3=3;
				if (c3=='T')
					c3=4;
				if (c3=='$')
					c3=5;
			}
			if (c4>5)
			{
				if (c4=='A')
					c4=1;
				if (c4=='C')
					c4=2;
				if (c4=='G')
					c4=3;
				if (c4=='T')
					c4=4;
				if (c4=='$')
					c4=5;
			}

			thedH += GetEnthalpy(c3,c4,c1,c2);
			thedS += GetEntropy(c3,c4,c1,c2);
		}

		return CalcTM(thedS,thedH);
	}
	inline float GetInitialEntropy() 
	{
		return 	-5.9f+rlogc;
	}

	static char getComplement(char,bool=false);
	static bool isMismatch(char,char);
private:
	float dH[6][6][6][6];  // A-C-G-T + gap + initiation (dangling end, $ sign)
	float dS[6][6][6][6];
};

#endif // !defined(AFX_NNPARAMS_H__05604705_84E8_11D4_A001_000000000000__INCLUDED_)
