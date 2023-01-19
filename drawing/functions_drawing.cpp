#include "functions_drawing.hpp"

float r = sqrt(2)-1;

float area2(float a, float b, float c){
    // area^2 of a triangle with edge lengths a b c
    return (a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c) / 16.;
}

float dihedral(float a, float b, float c, float d, float e, float f){
    // a=AD, b=BD, c=CD, d=BC, e=CA, f=AB
    //dihedral angle of edge a in the tetrahedron
    float Y2 = area2(a,c,e); //OCA
    float Z2 = area2(a,b,f); //OBA
    float H2 = (4.*pow(a,2)*pow(d,2) - pow((pow(b,2)+pow(e,2))-(pow(c,2)+pow(f,2)),2)) / 16.;
    return abs(acos( (-H2+Y2+Z2) / (2.*sqrt(Y2*Z2)) ));
}

float solid(float a, float b, float c, float d, float e, float f){
    // important: 3 first should share a vertex, the next ones are the opposites
    // a=AD, b=BD, c=CD, d=BC, e=CA, f=AB
    //solid angle of edges a,b,c in the tetrahedron (vertex D)
    return dihedral(a,b,c,d,e,f) + dihedral(b,d,f,e,a,c) + dihedral(c,d,e,f,a,b) - M_PI;    
}

float altitude(float ab, float ac, float ad, float bc, float bd, float cd){ // altitude of D 
    float xc=(pow(ab,2)+pow(ac,2)-pow(bc,2))/(float(2)*ab);
    float yc=sqrt(pow(ac,2)-pow(xc,2));
    float xd=(pow(ab,2)+pow(ad,2)-pow(bd,2))/(float(2)*ab);
    float yd=-(pow(ab,2) - pow(bd,2) + pow(cd,2) - pow(xc,2) - float(2)*ab*xd + float(2)*xc*xd - pow(yc,2))/(float(2)*yc);
    float zd=sqrt(pow(ad,2)-pow(xd,2)-pow(yd,2));
    return zd;
}


bool containsNoTriangle(float a,float b,float c){
    // true if the triple of intervals never satisfies triangle inequality
    return (a >= b+c || b >= a+c || c >= b+a);
}


// --------------------------------------- Tetr ---------------------------------------

string infos(float rA,float rB,float rC,float rD,float AB,float AC,float AD,float BC,float BD,float CD,bool detailed = false){
    string inf = "--- Tetr --- \n";
    inf += "  radii: " + to_string(rA) + " " + to_string(rB) +  " " + to_string(rC) +  " " + to_string(rD) + "\n";
    inf +=  "  edges: "+to_string(AB)+" "+to_string(AC)+" "+to_string(AD)+  " "+to_string(BC)+" "+to_string(BD)+" "+to_string(CD)+"\n";
    if (detailed){
	inf += "  volume: "+to_string(volume(AB,AC,AD,BC,BD,CD))+"\n";
	inf += "  cover: "+to_string(cover(rA,rB,rC,rD,AB,AC,AD,BC,BD,CD))+"\n";
    }
    inf += "  density: "+to_string(density(rA,rB,rC,rD,AB,AC,AD,BC,BD,CD))+"\n";
    if (detailed){
	inf += "  support sphere radius: "+to_string(supportR(rA,rB,rC,rD,AB,AC,AD,BC,BD,CD));
    }
    return inf;
}

float volume(float AB,float AC,float AD,float BC,float BD,float CD){
    float AD2 = pow(AD,2);
    float BD2 = pow(BD,2);
    float CD2 = pow(CD,2);
    float X = BD2 + CD2 - pow(BC,2);
    float Y = AD2 + CD2 - pow(AC,2);
    float Z = AD2 + BD2 - pow(AB,2);
    float V = sqrt(double(4)*AD2*BD2*CD2-AD2*pow(X,2)-BD2*pow(Y,2)-CD2*pow(Z,2)+X*Y*Z)/double(12);
    return V;
}


float cover(float rA,float rB,float rC,float rD,float AB,float AC,float AD,float BC,float BD,float CD){
    float covA = solid(AB, AC, AD, CD, BD, BC) * ((4./3.)*M_PI*pow(rA,3)) / (4.*M_PI);
    float covB = solid(AB, BC, BD, CD, AD, AC) * ((4./3.)*M_PI*pow(rB,3)) / (4.*M_PI);
    float covC = solid(AC, BC, CD, BD, AD, AB) * ((4./3.)*M_PI*pow(rC,3)) / (4.*M_PI);
    float covD = solid(AD, BD, CD, BC, AC, AB) * ((4./3.)*M_PI*pow(rD,3)) / (4.*M_PI);
    return covA+covB+covC+covD;
}


float density(float rA,float rB,float rC,float rD,float AB,float AC,float AD,float BC,float BD,float CD){
    return cover(rA,rB,rC,rD,AB,AC,AD,BC,BD,CD)/volume(AB,AC,AD,BC,BD,CD);
}


bool bad_radius(float ra, float rb, float rc, float rd, float ab, float ac, float ad, float bc, float bd, float cd){ // check if rx < zx for some x (sphere intersect a face)
    //    cout << leq(altitude(ab,ac,ad,bc,bd,cd),rd) <<  leq(altitude(bc,bd,ad,cd,ac,ad),ra) << leq(altitude(cd,ac,ab,ad,bd,ab),rb) << leq(altitude(ad,bd,cd,ab,ac,bc),rc) << endl;
    //cout << ra << " " << rb<< " " << rc<< " " << rd<< " " << ab<< " " << ac<< " " << ad<< " " << bc<< " " << bd<< " " << cd << endl;
    return (altitude(ab,ac,ad,bc,bd,cd) < rd ||  altitude(bc,bd,ab,cd,ac,ad) < ra || altitude(cd,ac,ab,ad,bd,ab) < rb || altitude(ad,bd,cd,ab,ac,bc) < rc);
}

float supportR(float ra, float rb, float rc, float rd, float ab, float ac, float ad, float bc, float bd, float cd){
    // attention, in certain cases if radius>1, we will say 1
    // coordinates of the vertices
    // A:(0,0,0), B:(ab,0,0), C:(xc,yc,0), D(xd,yd,zd)
    float xc=(pow(ab,2)+pow(ac,2)-pow(bc,2))/(float(2)*ab);
    float yc=sqrt(pow(ac,2)-pow(xc,2));
    float xd=(pow(ab,2)+pow(ad,2)-pow(bd,2))/(float(2)*ab);
    float yd=-(pow(ab,2) - pow(bd,2) + pow(cd,2) - pow(xc,2) - float(2)*ab*xd + float(2)*xc*xd - pow(yc,2))/(float(2)*yc);
    if (pow(ad,2)-pow(xd,2)-pow(yd,2) < 0) {
	return 100;
    }
    float zd=sqrt(pow(ad,2)-pow(xd,2)-pow(yd,2)); // 
    // the wanted radius is the  root of A*r^2+B*r+C=0
    float A=float(4)*(pow(ab,2)*pow(ra,2) - float(2)*pow(ab,2)*ra*rd + pow(ab,2)*pow(rd,2) + (pow(ra,2) - float(2)*ra*rb + pow(rb,2))*pow(xd,2) - float(2)*(ab*pow(ra,2) - ab*ra*rb - (ab*ra - ab*rb)*rd)*xd)*pow(yc,2) - float(8)*(pow(ab,2)*pow(ra,2) - pow(ab,2)*ra*rc - (pow(ab,2)*ra - pow(ab,2)*rc)*rd - (ab*pow(ra,2) - ab*ra*rb - (ab*ra - ab*rb)*rd)*xc - (ab*pow(ra,2) - ab*ra*rb - (ab*ra - ab*rb)*rc - (pow(ra,2) - float(2)*ra*rb + pow(rb,2))*xc)*xd)*yc*yd + float(4)*(pow(ab,2)*pow(ra,2) - float(2)*pow(ab,2)*ra*rc + pow(ab,2)*pow(rc,2) + (pow(ra,2) - float(2)*ra*rb + pow(rb,2))*pow(xc,2) - float(2)*(ab*pow(ra,2) - ab*ra*rb - (ab*ra - ab*rb)*rc)*xc)*pow(yd,2) + float(4)*(pow(ab,2)*pow(ra,2) - float(2)*pow(ab,2)*ra*rc + pow(ab,2)*pow(rc,2) + (pow(ra,2) - float(2)*ra*rb + pow(rb,2))*pow(xc,2) - (pow(ab,2) - pow(ra,2) + float(2)*ra*rb - pow(rb,2))*pow(yc,2) - float(2)*(ab*pow(ra,2) - ab*ra*rb - (ab*ra - ab*rb)*rc)*xc)*pow(zd,2);
    float B=-float(4)*(pow(ab,2)*ra - pow(ab,2)*rc - (ab*ra - ab*rb)*xc)*yc*pow(yd,3) + float(4)*(pow(ab,2)*pow(ra,3) - pow(ab,2)*pow(ra,2)*rd - pow(ab,2)*ra*pow(rd,2) + pow(ab,2)*pow(rd,3) - (ab*ra - ab*rb)*pow(xd,3) + (float(2)*pow(ab,2)*ra + pow(ra,3) - ra*pow(rb,2) + pow(rb,3) - pow(ab,2)*rd - (pow(ab,2) + pow(ra,2))*rb)*pow(xd,2) - (pow(ab,3)*ra + float(2)*ab*pow(ra,3) - ab*pow(ra,2)*rb - ab*ra*pow(rb,2) - (ab*ra - ab*rb)*pow(rd,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*rd)*xd)*pow(yc,2) + float(4)*(pow(ab,2)*pow(ra,3) - pow(ab,2)*pow(ra,2)*rc - pow(ab,2)*ra*pow(rc,2) + pow(ab,2)*pow(rc,3) - (ab*ra - ab*rb)*pow(xc,3) + (float(2)*pow(ab,2)*ra + pow(ra,3) - ra*pow(rb,2) + pow(rb,3) - pow(ab,2)*rc - (pow(ab,2) + pow(ra,2))*rb)*pow(xc,2) + (float(2)*pow(ab,2)*ra - pow(ab,2)*rc - pow(ab,2)*rd - (ab*ra - ab*rb)*xc - (ab*ra - ab*rb)*xd)*pow(yc,2) - (pow(ab,3)*ra + float(2)*ab*pow(ra,3) - ab*pow(ra,2)*rb - ab*ra*pow(rb,2) - (ab*ra - ab*rb)*pow(rc,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*rc)*xc)*pow(yd,2) + float(4)*(pow(ab,2)*pow(ra,3) - pow(ab,2)*pow(ra,2)*rc - pow(ab,2)*ra*pow(rc,2) + pow(ab,2)*pow(rc,3) - (ab*ra - ab*rb)*pow(xc,3) + (float(2)*pow(ab,2)*ra + pow(ra,3) - ra*pow(rb,2) + pow(rb,3) - pow(ab,2)*rc - (pow(ab,2) + pow(ra,2))*rb)*pow(xc,2) + (pow(ab,2)*ra + pow(ra,3) - ra*pow(rb,2) + pow(rb,3) - pow(ab,2)*rc - pow(ab,2)*rd - (pow(ab,2) + pow(ra,2))*rb - (ab*ra - ab*rb)*xc - (ab*ra - ab*rb)*xd)*pow(yc,2) - (pow(ab,2)*ra - pow(ab,2)*rc - (ab*ra - ab*rb)*xc)*yc*yd - (pow(ab,3)*ra + float(2)*ab*pow(ra,3) - ab*pow(ra,2)*rb - ab*ra*pow(rb,2) - (ab*ra - ab*rb)*pow(rc,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*rc)*xc)*pow(zd,2) - float(4)*((pow(ab,2)*ra - pow(ab,2)*rd - (ab*ra - ab*rb)*xd)*pow(yc,3) + (float(2)*pow(ab,2)*pow(ra,3) - pow(ab,2)*pow(ra,2)*rc - pow(ab,2)*ra*pow(rc,2) - (pow(ab,2)*ra - pow(ab,2)*rc)*pow(rd,2) + (pow(ab,2)*ra - pow(ab,2)*rd)*pow(xc,2) + (pow(ab,2)*ra - pow(ab,2)*rc - (ab*ra - ab*rb)*xc)*pow(xd,2) - (pow(ab,2)*pow(ra,2) - pow(ab,2)*pow(rc,2))*rd - (pow(ab,3)*ra + float(2)*ab*pow(ra,3) - ab*pow(ra,2)*rb - ab*ra*pow(rb,2) - (ab*ra - ab*rb)*pow(rd,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*rd)*xc - (pow(ab,3)*ra + float(2)*ab*pow(ra,3) - ab*pow(ra,2)*rb - ab*ra*pow(rb,2) - (ab*ra - ab*rb)*pow(rc,2) + (ab*ra - ab*rb)*pow(xc,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*rc - float(2)*(pow(ab,2)*ra + pow(ra,3) - ra*pow(rb,2) + pow(rb,3) - (pow(ab,2) + pow(ra,2))*rb)*xc)*xd)*yc)*yd;
    float C=pow(ab,2)*pow(yc,2)*pow(yd,4) + pow(ab,2)*pow(yc,2)*pow(zd,4) - float(2)*(pow(ab,2)*pow(yc,3) + (pow(ab,2)*pow(ra,2) - pow(ab,2)*pow(rc,2) + pow(ab,2)*pow(xc,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*xc)*yc)*pow(yd,3) + (pow(ab,2)*pow(ra,4) - float(2)*pow(ab,2)*pow(ra,2)*pow(rd,2) + pow(ab,2)*pow(rd,4) + pow(ab,2)*pow(xd,4) - float(2)*(pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*pow(xd,3) + (pow(ab,4) + float(4)*pow(ab,2)*pow(ra,2) + pow(ra,4) + pow(rb,4) - float(2)*pow(ab,2)*pow(rd,2) - float(2)*(pow(ab,2) + pow(ra,2))*pow(rb,2))*pow(xd,2) - float(2)*(pow(ab,3)*pow(ra,2) + ab*pow(ra,4) - ab*pow(ra,2)*pow(rb,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*pow(rd,2))*xd)*pow(yc,2) + (pow(ab,2)*pow(ra,4) - float(2)*pow(ab,2)*pow(ra,2)*pow(rc,2) + pow(ab,2)*pow(rc,4) + pow(ab,2)*pow(xc,4) + pow(ab,2)*pow(yc,4) - float(2)*(pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*pow(xc,3) + (pow(ab,4) + float(4)*pow(ab,2)*pow(ra,2) + pow(ra,4) + pow(rb,4) - float(2)*pow(ab,2)*pow(rc,2) - float(2)*(pow(ab,2) + pow(ra,2))*pow(rb,2))*pow(xc,2) + float(2)*(float(2)*pow(ab,2)*pow(ra,2) - pow(ab,2)*pow(rc,2) - pow(ab,2)*pow(rd,2) + pow(ab,2)*pow(xc,2) + pow(ab,2)*pow(xd,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*xc - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*xd)*pow(yc,2) - float(2)*(pow(ab,3)*pow(ra,2) + ab*pow(ra,4) - ab*pow(ra,2)*pow(rb,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*pow(rc,2))*xc)*pow(yd,2) + (pow(ab,2)*pow(ra,4) - float(2)*pow(ab,2)*pow(ra,2)*pow(rc,2) + pow(ab,2)*pow(rc,4) + pow(ab,2)*pow(xc,4) + pow(ab,2)*pow(yc,4) + float(2)*pow(ab,2)*pow(yc,2)*pow(yd,2) - float(2)*(pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*pow(xc,3) + (pow(ab,4) + float(4)*pow(ab,2)*pow(ra,2) + pow(ra,4) + pow(rb,4) - float(2)*pow(ab,2)*pow(rc,2) - float(2)*(pow(ab,2) + pow(ra,2))*pow(rb,2))*pow(xc,2) + (pow(ab,4) + float(2)*pow(ab,2)*pow(ra,2) + pow(ra,4) + pow(rb,4) - float(2)*pow(ab,2)*pow(rc,2) - float(2)*pow(ab,2)*pow(rd,2) + float(2)*pow(ab,2)*pow(xc,2) + float(2)*pow(ab,2)*pow(xd,2) - float(2)*(pow(ab,2) + pow(ra,2))*pow(rb,2) - float(2)*(pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*xc - float(2)*(pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*xd)*pow(yc,2) - float(2)*(pow(ab,3)*pow(ra,2) + ab*pow(ra,4) - ab*pow(ra,2)*pow(rb,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*pow(rc,2))*xc - float(2)*(pow(ab,2)*pow(yc,3) + (pow(ab,2)*pow(ra,2) - pow(ab,2)*pow(rc,2) + pow(ab,2)*pow(xc,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*xc)*yc)*yd)*pow(zd,2) - float(2)*((pow(ab,2)*pow(ra,2) - pow(ab,2)*pow(rd,2) + pow(ab,2)*pow(xd,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*xd)*pow(yc,3) + (pow(ab,2)*pow(ra,4) - pow(ab,2)*pow(ra,2)*pow(rc,2) - (pow(ab,2)*pow(ra,2) - pow(ab,2)*pow(rc,2))*pow(rd,2) + (pow(ab,2)*pow(ra,2) - pow(ab,2)*pow(rd,2))*pow(xc,2) + (pow(ab,2)*pow(ra,2) - pow(ab,2)*pow(rc,2) + pow(ab,2)*pow(xc,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*xc)*pow(xd,2) - (pow(ab,3)*pow(ra,2) + ab*pow(ra,4) - ab*pow(ra,2)*pow(rb,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*pow(rd,2))*xc - (pow(ab,3)*pow(ra,2) + ab*pow(ra,4) - ab*pow(ra,2)*pow(rb,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*pow(rc,2) + (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*pow(xc,2) - (pow(ab,4) + float(2)*pow(ab,2)*pow(ra,2) + pow(ra,4) + pow(rb,4) - float(2)*(pow(ab,2) + pow(ra,2))*pow(rb,2))*xc)*xd)*yc)*yd;
    float rs = 0;
    // usual formula when A<>0
    if (A*C != 0)
	{
	    float D=pow(B,2)-float(4)*A*C;
	    rs=(C<0 ? (-B+sqrt(D))/(2*A) : (-B-sqrt(D))/(2*A));
	}
    // alternative formula (necessary if A contains 0, useful if A very small)
    float x=4*A*C/pow(B,2);
    if (x<0.78){
	rs=-C/B*(1+x);
    }

    // if A,B contain 0 ? 
    //if (-C / (A+B) > 1){
	
    //}
    
    return rs;
} 

bool containsNoTetr(float AB,float AC,float AD,float BC,float BD,float CD){ // not sure it is enough
    bool noValidFaces = containsNoTriangle(AB,BC,AC) || containsNoTriangle(AB,BD,AD) || containsNoTriangle(AC,CD,AD)  || containsNoTriangle(BC,CD,BD);
    return noValidFaces;
    /* if (noValidFaces){return true;}
    else {return false;} // forget about the volume
    float AD2 = pow(AD,2);
    float BD2 = pow(BD,2);
    float CD2 = pow(CD,2);
    float X = BD2 + CD2 - pow(BC,2);
    float Y = AD2 + CD2 - pow(AC,2);
    float Z = AD2 + BD2 - pow(AB,2);
    float V = double(4)*AD2*BD2*CD2-AD2*pow(X,2)-BD2*pow(Y,2)-CD2*pow(Z,2)+X*Y*Z;
    bool noVol = (upper(V) <= 0);
    return noVol;    */
}

bool containsNoFM(float rA,float rB,float rC,float rD,float AB,float AC,float AD,float BC,float BD,float CD){
    // True if there are no FM tetrahedra in the set of this t
    return (containsNoTetr(AB,AC,AD,BC,BD,CD) || bad_radius(rA,rB,rC,rD,AB,AC,AD,BC,BD,CD) || supportR(rA,rB,rC,rD,AB,AC,AD,BC,BD,CD) > r);
}
      
