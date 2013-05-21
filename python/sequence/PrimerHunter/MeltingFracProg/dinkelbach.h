//Module: Achievement of the Dinkelbach algorithm

#ifndef __dinkelbach_h
#define __dinkelbach_h

#include <iostream>
#include <strstream>
#include <cstring>
#include "nnparams.h"
#include "galign.h"

using namespace std;

//typedef class dinkelbach* DINKEL;

class dinkelbach
{
  private:
    PNNParams myDinkParams;
    GGAlign myGAlign;

  public:
    dinkelbach(PNNParams, GGAlign);
    ~dinkelbach();
    int iteration(float);
};

#endif /* __dinkelbach_h */
