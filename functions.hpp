#ifndef TETRHPP_INCLUDED
#define TETRHPP_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/io.hpp>
#include <string.h>
#include <math.h>

using namespace boost::numeric;
using namespace interval_lib;
using namespace std;

typedef interval<double, policies<save_state<rounded_transc_std<double>>, checking_base<double>>> I;

extern I PI;
extern I r;

string I2str(I x);

bool leq(I a, I b);

string infos(I rA,I rB,I rC,I rD,I AB,I AC,I AD,I BC,I BD,I CD);
I volume(I AB,I AC,I AD,I BC,I BD,I CD);
I cover(I rA,I rB,I rC,I rD,I AB,I AC,I AD,I BC,I BD,I CD);
I density(I rA,I rB,I rC,I rD,I AB,I AC,I AD,I BC,I BD,I CD);
bool bad_radius(I rA,I rB,I rC,I rD,I AB,I AC,I AD,I BC,I BD,I CD);
I supportR(I rA,I rB,I rC,I rD,I AB,I AC,I AD,I BC,I BD,I CD);
bool containsNoTetr(I AB,I AC,I AD,I BC,I BD,I CD); // if edge inequalities are not satisfied
bool containsNoFM(I rA,I rB,I rC,I rD,I AB,I AC,I AD,I BC,I BD,I CD); // if not FM


#endif
