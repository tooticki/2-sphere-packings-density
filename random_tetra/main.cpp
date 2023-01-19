#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <bits/stdc++.h> 
#include <boost/numeric/interval.hpp> // Library for interval arithmetic
#include <boost/numeric/interval/io.hpp>
#include <deque>
#include "functions.hpp" // useful functions

using namespace boost::numeric;
using namespace interval_lib;
using namespace std;
using namespace std::chrono;


typedef interval<double, policies<save_state<rounded_transc_std<double>>, checking_base<double>>> I;

I PI = pi<I>();
I r = I(sqrt(2)-1);
I R = I(1);
double eEPS = 0.000000000001;
I AB = 0;
I RA, RB, RC, RD; 

class bad_tetr_exception: public exception{
  virtual const char* what() const throw()  {
    return "Found a bad tetrahedron";
  }
} bad_tetr;

struct Block{
    I edges[5];
};

int N_blocks = 0;
int n_del_blocks = 0;
int n_del[2] = {0,0};
deque <Block> blocks;
int BFSlimit = 4194304; //1048576 for magi
int DFSlimit = BFSlimit/2;


//division in two of the greatest dimension

//with spheres of radii 1,0.414213562373095?,1,1
//with edges of lengths 1.41421356237310?, 2.00000000000000?, 2.00000000000000?, 1.4142135623731?, 1.41421356237310?, 2.69245416472622?

//stretched 11r1: Vol=0.27492393477375?, E=-0.00534375425043319, Density=0.812542027810835
//Tetr* t_tightest = new Tetr(?) 
I ud = I(0.812542027810795, 0.812542027810875); // upper_density (enough precise for now)
double lower_ud = 0.812542027810795, upper_ud = 0.812542027810875;

// ------------------------- Generalization ------------------------------


/* possible symmetries:
 1: rA==rB and rC==rD => AC<=AD,AC<=BC,AC<=BD 
 2: rA!=rB and rC==rD => AC<=AD 
 3: rA==rB and rC!=rD => AC<=BC
*/

void print_block(Block b){
    for(int i = 0; i<5; i++){
	cout << b.edges[i] << " ";
    }
    cout << endl;
}

int run_blocks(I ac, I ad, I bc, I bd, I cd){
    Block b;
    I e[5] = {ac, ad, bc, bd, cd};
    memcpy(b.edges, e, 5*sizeof(I));    
    blocks.push_back(b);
    
    I AC, AD, BC, BD, CD;
    I dens;
    int maxi = 0;
    double maxw = 0;	    
    Block bb0;
    Block bb1;
    bool bf = true; // bfs if true otherwise dfs(to limit memory)
    
    while(blocks.size()>0){
	if (bf){	    
	    if (double(blocks.size()/max(N_blocks,1)) >= 2 || double(N_blocks/(blocks.size())) >= 2 ){
		N_blocks = blocks.size();
		cout << N_blocks << " blocks;  " << n_del[0]+n_del[1]+n_del[2] << " deleted; width: " << maxw << endl;
	    }
	    if(blocks.size() >= BFSlimit){
		bf = false;            
		cout << "Passing to DFS" << endl;
	    }
	}
	else{
	    if(int(blocks.size())-N_blocks >= DFSlimit/32|| int(blocks.size())-N_blocks <= -DFSlimit/32){
		//return 0 ; // to comment
		N_blocks = blocks.size();
		cout << N_blocks << " blocks;  " << n_del[0]+n_del[1]+n_del[2] << " deleted; width: " << maxw  << endl;		
	    }
	    if (blocks.size() <= DFSlimit ){
		bf = true;
		cout << "Repassing to BFS" << endl;
	    }
	}
	b = blocks.front();
	blocks.pop_front();

	AC=b.edges[0]; AD=b.edges[1]; BC=b.edges[2]; BD=b.edges[3]; CD=b.edges[4];
	dens = density(RA, RB, RC, RD, AB, AC, AD ,BC, BD, CD);
	
	if(upper(dens)<lower_ud){// good	    
	    n_del[0]++;
	}
	else if(lower(AC)>upper(BC) || lower(AC)>upper(AD) || lower(AC)>upper(BD)){ // AC = min (BC, AD, BD), sym only for 1111 and 11rr and rr11
	    n_del[1]++;
	}
	else if(containsNoTetr(AB,AC,AD,BC,BD,CD) || bad_radius(RA,RB,RC,RD,AB,AC,AD,BC,BD,CD) || lower(supportR(RA,RB,RC,RD,AB,AC,AD,BC,BD,CD)) > upper(r)){ // no FM-tetr
	    n_del[1]++;
	}
	else{ // subdivide and add to queue
	    maxi = 0;
	    maxw = 0.;
	    for(int i = 0; i<5; i++){
		bb0.edges[i] = b.edges[i];
		bb1.edges[i] = b.edges[i];
		if (width(b.edges[i])>maxw){
		    maxi = i;
		    maxw = width(b.edges[i]);		    
		}
	    }
	    bb0.edges[maxi] = bisect(b.edges[maxi]).first;
	    bb1.edges[maxi] = bisect(b.edges[maxi]).second;
	    if(bf){
		blocks.push_back(bb0);
		blocks.push_back(bb1);
	    }
	    else{		
		blocks.push_front(bb0);
		blocks.push_front(bb1);
	    }
	}
	if(upper(dens) < 1 && lower(dens) > upper_ud){ // BAD
	    throw bad_tetr;
	}
    }
    return 0;
}



I edge_length(I r1,I r2, double EPS){
    return I(lower(r1+r2+EPS), upper(r1+r2+r*2.));
}


int test_tetrahedra(I rA, I rB, I rC, I rD, int c_type){
    // if c_type = 1:  just a contact between rA and rB (AB=rA+rB)
    // if c_type = 21:  additional contact rA rC (AC=rA+rC)
    // if c_type = 22:  additional contact rA rD (AD=rA+rD)
    // if c_type = 23:  additional contact rB rC (BC=rB+rC)
    // if c_type = 24:  additional contact rB rD (BD=rB+rD)
    // if c_type = 25:  additional contact rC rD (CD=rC+rD)
    double EPS = 0.; // we ignore EPS-tight for now but 0 for 1111
    AB = rA+rB;
    I AC=edge_length(rA,rC,EPS), AD=edge_length(rA,rD,EPS), BC=edge_length(rB,rC,EPS);
    I BD=edge_length(rB,rD,EPS), CD=edge_length(rC,rD,EPS);
    
    switch(c_type){
    case 21:
	AC = rA+rC; break;	
    case 22:
	AD = rA+rD; break;
    case 23:
	BC = rB+rC; break;
    case 24:
	BD = rB+rD; break;
    case 25:
	CD = rC+rD; break;
    }

    RA=rA;
    RB=rB;
    RC=rC;
    RD=rD;
    
    cout << "testing all tetrahedra in "<<AB<<" "<<AC<<" "<<AD<<" "<<BC<<" "<<BD<<" "<<CD<<endl;
    run_blocks(AC,AD,BC,BD,CD);
    cout << "DONE "<< endl;
    return 1;

}

I rval (char* c){
    return ((c[0] == '1') ? R : r);
}

int main(int argc, char *argv[])
{
    cout << "upper density: " <<  density(R,r,R,R,R+r,2.*R, 2.*R, R+r, R+r, 2.69245416472622) << endl;
    if (argc != 6){
	cout << "> ./program r1 r2 r3 r4 \n";
	return 0;
    }
    cout << "Radii : ";
    for (int i = 1; i < 5; i++){
	cout << rval(argv[i]) << " ";
    }
    //contacts: 12, 13, 14
    I rA = argv[1], rB = argv[2], rC = argv[3], rD = argv[4];
    I AB = rA+rB, AC=rA+rC, AD=rA+rD;
    I BC, BD, CD;
    for (i=0; i<50; i++){
	BC = rB+rC+//TODO
    }
    
    
} 


    //rrrr done with 1111 : 
    
    //11rr
    //1r1r
    //rr11

    //1rrr
    //rr1r
    
    //111r (hardest)
    //1r11
    // TODO (EPS!=0)



    /*
    Tetr* t = new Tetr(1,1,1,1,2,2,2*sqrt(2),2*sqrt(2),2,2);
    
    density
    cout << "DENS: " << t->infos(true) << t->containsNoTetr() << t->containsNoFM() << endl;*/
    //cout << density(1,1,1,1,2,2,2*sqrt(2),2*sqrt(2),2,2);
