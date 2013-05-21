#include <iostream>
#include <strstream>
#include <cstring>
#include "nnparams.h"
#include "galign.h"
#include "dinkelbach.h"

using namespace std;


dinkelbach::dinkelbach(PNNParams myDiParams, GGAlign mGAlign)
{
  myDinkParams = myDiParams;
  myGAlign = mGAlign;
}


dinkelbach::~dinkelbach()
{
  //does nothing!
}

int dinkelbach::iteration(float TempK)
{
  int iteration = 0;
  myDinkParams->AlterTM(TempK);
  float precursor_oldTM = TempK;
  float newTM = TempK;
  float oldTM = TempK, lambda, oldlambda;
  lambda = 999999999999999999.9;
  int positive=0;

  do
  {
    iteration++;
    myGAlign->InitBorder();
    myGAlign->CalculateTable();
#ifdef _output_alignment
    myGAlign->OutputLocalAlignment(std::cout);
    cout << "Position i: " << myGAlign->maxloci << " Position j: " << myGAlign->maxlocj << " Melting temp: "<< newTM<<endl;
#endif
    oldlambda=lambda;
    lambda =myGAlign->GetFreeEnergyK(myGAlign->maxloci,myGAlign->maxlocj,newTM);
#ifdef _output_alignment
    cout << "Iteration " << iteration << " lambda " << lambda << endl;
#endif
    if(lambda < 0.0) 
    {
	if(!positive) 
	{
		//Initial melting temperature greater than expected
		newTM = myGAlign->GetEnthalpy(myGAlign->maxloci,myGAlign->maxlocj)/myGAlign->GetEntropy(myGAlign->maxloci,myGAlign->maxlocj)-1;
		lambda = oldlambda;
    		myDinkParams->AlterTM(newTM);
	} 
	else 
	{
		precursor_oldTM = oldTM;
	}
    } 
    else 
    {
	positive = 1;
    	oldTM = newTM;
    	newTM = myGAlign->GetMeltingTempK(myGAlign->maxloci,myGAlign->maxlocj);
    	myDinkParams->AlterTM(newTM);
    }
  }
  while(!positive || (newTM - oldTM > 0.001  && (lambda>0.0)&&(lambda<oldlambda)));

  if(newTM - oldTM > 0.001 && ((lambda<0.0)||(lambda>oldlambda)))
  {
    iteration++;
    myDinkParams->AlterTM(precursor_oldTM);
    myGAlign->InitBorder();
    myGAlign->CalculateTable();
#ifdef _output_alignment
    myGAlign->OutputLocalAlignment(std::cout);
#endif
    lambda = myGAlign->GetFreeEnergyK(myGAlign->maxloci,myGAlign->maxlocj,precursor_oldTM);
    newTM = myGAlign->GetMeltingTempK(myGAlign->maxloci,myGAlign->maxlocj);
    myDinkParams->AlterTM(newTM);
#ifdef _output_alignment
    cout << "Iteration " << iteration << " lambda " << lambda << " extra if" << endl;
#endif
  }
#ifdef _output_alignment
  cout << "new value for TM in K: " << newTM << endl;
  cout << "new value for TM in °C: " << newTM - 273.15 << endl;
#endif
  return iteration;
}
