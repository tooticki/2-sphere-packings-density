#ifndef TETRHPP_floatNCLUDED
#define TETRHPP_floatNCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <math.h>

using namespace std;


string float2str(float x);

bool leq(float a, float b);

string infos(float rA,float rB,float rC,float rD,float AB,float AC,float AD,float BC,float BD,float CD);
float volume(float AB,float AC,float AD,float BC,float BD,float CD);
float cover(float rA,float rB,float rC,float rD,float AB,float AC,float AD,float BC,float BD,float CD);
float density(float rA,float rB,float rC,float rD,float AB,float AC,float AD,float BC,float BD,float CD);
bool bad_radius(float rA,float rB,float rC,float rD,float AB,float AC,float AD,float BC,float BD,float CD);
float supportR(float rA,float rB,float rC,float rD,float AB,float AC,float AD,float BC,float BD,float CD);
bool containsNoTetr(float AB,float AC,float AD,float BC,float BD,float CD); // if edge inequalities are not satisfied
bool containsNoFM(float rA,float rB,float rC,float rD,float AB,float AC,float AD,float BC,float BD,float CD); // if not FM


#endif
