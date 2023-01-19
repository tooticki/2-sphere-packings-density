#include "functions.hpp"


bool leq (I a, I b){ // a<b
    return upper(a) < lower(b);
}

string I2str(I x){
    return '[' + to_string(x.lower()) + ',' + to_string(x.upper()) + ']';
}

I area2(I a, I b, I c){
    // area^2 of a triangle with edge lengths a b c
    return (a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c) / 16.;
}

I dihedral(I a, I b, I c, I d, I e, I f){
    // a=AD, b=BD, c=CD, d=BC, e=CA, f=AB
    //dihedral angle of edge a in the tetrahedron
    I Y2 = area2(a,c,e); //OCA
    I Z2 = area2(a,b,f); //OBA
    I H2 = (4.*pow(a,2)*pow(d,2) - pow((pow(b,2)+pow(e,2))-(pow(c,2)+pow(f,2)),2)) / 16.;
    return abs(acos( (-H2+Y2+Z2) / (2.*sqrt(Y2*Z2)) ));
}

I solid(I a, I b, I c, I d, I e, I f){
    // important: 3 first should share a vertex, the next ones are the opposites
    // a=AD, b=BD, c=CD, d=BC, e=CA, f=AB
    //solid angle of edges a,b,c in the tetrahedron (vertex D)
    return dihedral(a,b,c,d,e,f) + dihedral(b,d,f,e,a,c) + dihedral(c,d,e,f,a,b) - PI;    
}

I altitude(I ab, I ac, I ad, I bc, I bd, I cd){ // altitude of D 
    I xc=(pow(ab,2)+pow(ac,2)-pow(bc,2))/(I(2)*ab);
    I yc=sqrt(pow(ac,2)-pow(xc,2));
    I xd=(pow(ab,2)+pow(ad,2)-pow(bd,2))/(I(2)*ab);
    I yd=-(pow(ab,2) - pow(bd,2) + pow(cd,2) - pow(xc,2) - I(2)*ab*xd + I(2)*xc*xd - pow(yc,2))/(I(2)*yc);
    I zd=sqrt(pow(ad,2)-pow(xd,2)-pow(yd,2));
    return zd;
}


bool containsNoTriangle(I a,I b,I c){
    // true if the triple of intervals never satisfies triangle inequality
    return (lower(a) > upper(b+c) || lower(b) > upper(a+c) || lower(c) > upper(b+a));
}


// --------------------------------------- Tetr ---------------------------------------

string infos(I rA,I rB,I rC,I rD,I AB,I AC,I AD,I BC,I BD,I CD,bool detailed = false){
    string inf = "--- Tetr --- \n";
    inf += "  radii: " + I2str(rA) + " " + I2str(rB) +  " " + I2str(rC) +  " " + I2str(rD) + "\n";
    inf +=  "  edges: "+I2str(AB)+" "+I2str(AC)+" "+I2str(AD)+  " "+I2str(BC)+" "+I2str(BD)+" "+I2str(CD)+"\n";
    if (detailed){
	inf += "  volume: "+I2str(volume(AB,AC,AD,BC,BD,CD))+"\n";
	inf += "  cover: "+I2str(cover(rA,rB,rC,rD,AB,AC,AD,BC,BD,CD))+"\n";
    }
    inf += "  density: "+I2str(density(rA,rB,rC,rD,AB,AC,AD,BC,BD,CD))+"\n";
    if (detailed){
	inf += "  support sphere radius: "+I2str(supportR(rA,rB,rC,rD,AB,AC,AD,BC,BD,CD));
    }
    return inf;
}

I volume(I AB,I AC,I AD,I BC,I BD,I CD){
    I AD2 = pow(AD,2);
    I BD2 = pow(BD,2);
    I CD2 = pow(CD,2);
    I X = BD2 + CD2 - pow(BC,2);
    I Y = AD2 + CD2 - pow(AC,2);
    I Z = AD2 + BD2 - pow(AB,2);
    I V = sqrt(double(4)*AD2*BD2*CD2-AD2*pow(X,2)-BD2*pow(Y,2)-CD2*pow(Z,2)+X*Y*Z)/double(12);
    return V;
}


I cover(I rA,I rB,I rC,I rD,I AB,I AC,I AD,I BC,I BD,I CD){
    I covA = solid(AB, AC, AD, CD, BD, BC) * ((4./3.)*PI*pow(rA,3)) / (4.*PI);
    I covB = solid(AB, BC, BD, CD, AD, AC) * ((4./3.)*PI*pow(rB,3)) / (4.*PI);
    I covC = solid(AC, BC, CD, BD, AD, AB) * ((4./3.)*PI*pow(rC,3)) / (4.*PI);
    I covD = solid(AD, BD, CD, BC, AC, AB) * ((4./3.)*PI*pow(rD,3)) / (4.*PI);
    return covA+covB+covC+covD;
}


I density(I rA,I rB,I rC,I rD,I AB,I AC,I AD,I BC,I BD,I CD){
    return cover(rA,rB,rC,rD,AB,AC,AD,BC,BD,CD)/volume(AB,AC,AD,BC,BD,CD);
}


bool bad_radius(I ra, I rb, I rc, I rd, I ab, I ac, I ad, I bc, I bd, I cd){ // check if rx < zx for some x (sphere intersect a face)
    //    cout << leq(altitude(ab,ac,ad,bc,bd,cd),rd) <<  leq(altitude(bc,bd,ad,cd,ac,ad),ra) << leq(altitude(cd,ac,ab,ad,bd,ab),rb) << leq(altitude(ad,bd,cd,ab,ac,bc),rc) << endl;
    //cout << ra << " " << rb<< " " << rc<< " " << rd<< " " << ab<< " " << ac<< " " << ad<< " " << bc<< " " << bd<< " " << cd << endl;
    return leq(altitude(ab,ac,ad,bc,bd,cd),rd) ||  leq(altitude(bc,bd,ab,cd,ac,ad),ra) || leq(altitude(cd,ac,ab,ad,bd,ab),rb) || leq(altitude(ad,bd,cd,ab,ac,bc),rc);
}

I supportR(I ra, I rb, I rc, I rd, I ab, I ac, I ad, I bc, I bd, I cd){
    // attention, in certain cases if radius>1, we will say 1
    // coordinates of the vertices
    // A:(0,0,0), B:(ab,0,0), C:(xc,yc,0), D(xd,yd,zd)
    I xc=(pow(ab,2)+pow(ac,2)-pow(bc,2))/(I(2)*ab);
    I yc=sqrt(pow(ac,2)-pow(xc,2));
    I xd=(pow(ab,2)+pow(ad,2)-pow(bd,2))/(I(2)*ab);
    I yd=-(pow(ab,2) - pow(bd,2) + pow(cd,2) - pow(xc,2) - I(2)*ab*xd + I(2)*xc*xd - pow(yc,2))/(I(2)*yc);
    I zd=sqrt(pow(ad,2)-pow(xd,2)-pow(yd,2)); // 
    // the wanted radius is the  root of A*r^2+B*r+C=0
    I A=I(4)*(pow(ab,2)*pow(ra,2) - I(2)*pow(ab,2)*ra*rd + pow(ab,2)*pow(rd,2) + (pow(ra,2) - I(2)*ra*rb + pow(rb,2))*pow(xd,2) - I(2)*(ab*pow(ra,2) - ab*ra*rb - (ab*ra - ab*rb)*rd)*xd)*pow(yc,2) - I(8)*(pow(ab,2)*pow(ra,2) - pow(ab,2)*ra*rc - (pow(ab,2)*ra - pow(ab,2)*rc)*rd - (ab*pow(ra,2) - ab*ra*rb - (ab*ra - ab*rb)*rd)*xc - (ab*pow(ra,2) - ab*ra*rb - (ab*ra - ab*rb)*rc - (pow(ra,2) - I(2)*ra*rb + pow(rb,2))*xc)*xd)*yc*yd + I(4)*(pow(ab,2)*pow(ra,2) - I(2)*pow(ab,2)*ra*rc + pow(ab,2)*pow(rc,2) + (pow(ra,2) - I(2)*ra*rb + pow(rb,2))*pow(xc,2) - I(2)*(ab*pow(ra,2) - ab*ra*rb - (ab*ra - ab*rb)*rc)*xc)*pow(yd,2) + I(4)*(pow(ab,2)*pow(ra,2) - I(2)*pow(ab,2)*ra*rc + pow(ab,2)*pow(rc,2) + (pow(ra,2) - I(2)*ra*rb + pow(rb,2))*pow(xc,2) - (pow(ab,2) - pow(ra,2) + I(2)*ra*rb - pow(rb,2))*pow(yc,2) - I(2)*(ab*pow(ra,2) - ab*ra*rb - (ab*ra - ab*rb)*rc)*xc)*pow(zd,2);
    I B=-I(4)*(pow(ab,2)*ra - pow(ab,2)*rc - (ab*ra - ab*rb)*xc)*yc*pow(yd,3) + I(4)*(pow(ab,2)*pow(ra,3) - pow(ab,2)*pow(ra,2)*rd - pow(ab,2)*ra*pow(rd,2) + pow(ab,2)*pow(rd,3) - (ab*ra - ab*rb)*pow(xd,3) + (I(2)*pow(ab,2)*ra + pow(ra,3) - ra*pow(rb,2) + pow(rb,3) - pow(ab,2)*rd - (pow(ab,2) + pow(ra,2))*rb)*pow(xd,2) - (pow(ab,3)*ra + I(2)*ab*pow(ra,3) - ab*pow(ra,2)*rb - ab*ra*pow(rb,2) - (ab*ra - ab*rb)*pow(rd,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*rd)*xd)*pow(yc,2) + I(4)*(pow(ab,2)*pow(ra,3) - pow(ab,2)*pow(ra,2)*rc - pow(ab,2)*ra*pow(rc,2) + pow(ab,2)*pow(rc,3) - (ab*ra - ab*rb)*pow(xc,3) + (I(2)*pow(ab,2)*ra + pow(ra,3) - ra*pow(rb,2) + pow(rb,3) - pow(ab,2)*rc - (pow(ab,2) + pow(ra,2))*rb)*pow(xc,2) + (I(2)*pow(ab,2)*ra - pow(ab,2)*rc - pow(ab,2)*rd - (ab*ra - ab*rb)*xc - (ab*ra - ab*rb)*xd)*pow(yc,2) - (pow(ab,3)*ra + I(2)*ab*pow(ra,3) - ab*pow(ra,2)*rb - ab*ra*pow(rb,2) - (ab*ra - ab*rb)*pow(rc,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*rc)*xc)*pow(yd,2) + I(4)*(pow(ab,2)*pow(ra,3) - pow(ab,2)*pow(ra,2)*rc - pow(ab,2)*ra*pow(rc,2) + pow(ab,2)*pow(rc,3) - (ab*ra - ab*rb)*pow(xc,3) + (I(2)*pow(ab,2)*ra + pow(ra,3) - ra*pow(rb,2) + pow(rb,3) - pow(ab,2)*rc - (pow(ab,2) + pow(ra,2))*rb)*pow(xc,2) + (pow(ab,2)*ra + pow(ra,3) - ra*pow(rb,2) + pow(rb,3) - pow(ab,2)*rc - pow(ab,2)*rd - (pow(ab,2) + pow(ra,2))*rb - (ab*ra - ab*rb)*xc - (ab*ra - ab*rb)*xd)*pow(yc,2) - (pow(ab,2)*ra - pow(ab,2)*rc - (ab*ra - ab*rb)*xc)*yc*yd - (pow(ab,3)*ra + I(2)*ab*pow(ra,3) - ab*pow(ra,2)*rb - ab*ra*pow(rb,2) - (ab*ra - ab*rb)*pow(rc,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*rc)*xc)*pow(zd,2) - I(4)*((pow(ab,2)*ra - pow(ab,2)*rd - (ab*ra - ab*rb)*xd)*pow(yc,3) + (I(2)*pow(ab,2)*pow(ra,3) - pow(ab,2)*pow(ra,2)*rc - pow(ab,2)*ra*pow(rc,2) - (pow(ab,2)*ra - pow(ab,2)*rc)*pow(rd,2) + (pow(ab,2)*ra - pow(ab,2)*rd)*pow(xc,2) + (pow(ab,2)*ra - pow(ab,2)*rc - (ab*ra - ab*rb)*xc)*pow(xd,2) - (pow(ab,2)*pow(ra,2) - pow(ab,2)*pow(rc,2))*rd - (pow(ab,3)*ra + I(2)*ab*pow(ra,3) - ab*pow(ra,2)*rb - ab*ra*pow(rb,2) - (ab*ra - ab*rb)*pow(rd,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*rd)*xc - (pow(ab,3)*ra + I(2)*ab*pow(ra,3) - ab*pow(ra,2)*rb - ab*ra*pow(rb,2) - (ab*ra - ab*rb)*pow(rc,2) + (ab*ra - ab*rb)*pow(xc,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*rc - I(2)*(pow(ab,2)*ra + pow(ra,3) - ra*pow(rb,2) + pow(rb,3) - (pow(ab,2) + pow(ra,2))*rb)*xc)*xd)*yc)*yd;
    I C=pow(ab,2)*pow(yc,2)*pow(yd,4) + pow(ab,2)*pow(yc,2)*pow(zd,4) - I(2)*(pow(ab,2)*pow(yc,3) + (pow(ab,2)*pow(ra,2) - pow(ab,2)*pow(rc,2) + pow(ab,2)*pow(xc,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*xc)*yc)*pow(yd,3) + (pow(ab,2)*pow(ra,4) - I(2)*pow(ab,2)*pow(ra,2)*pow(rd,2) + pow(ab,2)*pow(rd,4) + pow(ab,2)*pow(xd,4) - I(2)*(pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*pow(xd,3) + (pow(ab,4) + I(4)*pow(ab,2)*pow(ra,2) + pow(ra,4) + pow(rb,4) - I(2)*pow(ab,2)*pow(rd,2) - I(2)*(pow(ab,2) + pow(ra,2))*pow(rb,2))*pow(xd,2) - I(2)*(pow(ab,3)*pow(ra,2) + ab*pow(ra,4) - ab*pow(ra,2)*pow(rb,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*pow(rd,2))*xd)*pow(yc,2) + (pow(ab,2)*pow(ra,4) - I(2)*pow(ab,2)*pow(ra,2)*pow(rc,2) + pow(ab,2)*pow(rc,4) + pow(ab,2)*pow(xc,4) + pow(ab,2)*pow(yc,4) - I(2)*(pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*pow(xc,3) + (pow(ab,4) + I(4)*pow(ab,2)*pow(ra,2) + pow(ra,4) + pow(rb,4) - I(2)*pow(ab,2)*pow(rc,2) - I(2)*(pow(ab,2) + pow(ra,2))*pow(rb,2))*pow(xc,2) + I(2)*(I(2)*pow(ab,2)*pow(ra,2) - pow(ab,2)*pow(rc,2) - pow(ab,2)*pow(rd,2) + pow(ab,2)*pow(xc,2) + pow(ab,2)*pow(xd,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*xc - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*xd)*pow(yc,2) - I(2)*(pow(ab,3)*pow(ra,2) + ab*pow(ra,4) - ab*pow(ra,2)*pow(rb,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*pow(rc,2))*xc)*pow(yd,2) + (pow(ab,2)*pow(ra,4) - I(2)*pow(ab,2)*pow(ra,2)*pow(rc,2) + pow(ab,2)*pow(rc,4) + pow(ab,2)*pow(xc,4) + pow(ab,2)*pow(yc,4) + I(2)*pow(ab,2)*pow(yc,2)*pow(yd,2) - I(2)*(pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*pow(xc,3) + (pow(ab,4) + I(4)*pow(ab,2)*pow(ra,2) + pow(ra,4) + pow(rb,4) - I(2)*pow(ab,2)*pow(rc,2) - I(2)*(pow(ab,2) + pow(ra,2))*pow(rb,2))*pow(xc,2) + (pow(ab,4) + I(2)*pow(ab,2)*pow(ra,2) + pow(ra,4) + pow(rb,4) - I(2)*pow(ab,2)*pow(rc,2) - I(2)*pow(ab,2)*pow(rd,2) + I(2)*pow(ab,2)*pow(xc,2) + I(2)*pow(ab,2)*pow(xd,2) - I(2)*(pow(ab,2) + pow(ra,2))*pow(rb,2) - I(2)*(pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*xc - I(2)*(pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*xd)*pow(yc,2) - I(2)*(pow(ab,3)*pow(ra,2) + ab*pow(ra,4) - ab*pow(ra,2)*pow(rb,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*pow(rc,2))*xc - I(2)*(pow(ab,2)*pow(yc,3) + (pow(ab,2)*pow(ra,2) - pow(ab,2)*pow(rc,2) + pow(ab,2)*pow(xc,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*xc)*yc)*yd)*pow(zd,2) - I(2)*((pow(ab,2)*pow(ra,2) - pow(ab,2)*pow(rd,2) + pow(ab,2)*pow(xd,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*xd)*pow(yc,3) + (pow(ab,2)*pow(ra,4) - pow(ab,2)*pow(ra,2)*pow(rc,2) - (pow(ab,2)*pow(ra,2) - pow(ab,2)*pow(rc,2))*pow(rd,2) + (pow(ab,2)*pow(ra,2) - pow(ab,2)*pow(rd,2))*pow(xc,2) + (pow(ab,2)*pow(ra,2) - pow(ab,2)*pow(rc,2) + pow(ab,2)*pow(xc,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*xc)*pow(xd,2) - (pow(ab,3)*pow(ra,2) + ab*pow(ra,4) - ab*pow(ra,2)*pow(rb,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*pow(rd,2))*xc - (pow(ab,3)*pow(ra,2) + ab*pow(ra,4) - ab*pow(ra,2)*pow(rb,2) - (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*pow(rc,2) + (pow(ab,3) + ab*pow(ra,2) - ab*pow(rb,2))*pow(xc,2) - (pow(ab,4) + I(2)*pow(ab,2)*pow(ra,2) + pow(ra,4) + pow(rb,4) - I(2)*(pow(ab,2) + pow(ra,2))*pow(rb,2))*xc)*xd)*yc)*yd;
    I rs=I(-INFINITY,+INFINITY);
    // usual formula when A<>0
    if (not zero_in(A*C))
    {
        I D=square(B)-I(4)*A*C;
        rs=intersect(rs,(upper(C)<0 ? (-B+sqrt(D))/(I(2)*A) : (-B-sqrt(D))/(I(2)*A)));
    }
    // alternative formula (necessary if A contains 0, useful if A very small)
    I x=I(4)*A*C/square(B);
    x=hull(x,I(0)); // the interval x must contain 0
    if (upper(x)<0.78) rs=intersect(rs,-C/B*(I(1)+x));

    // if A,B contain 0 ? 
    //if (-C / (A+B) > 1){
	
    //}
    
    return rs;
} 

bool containsNoTetr(I AB,I AC,I AD,I BC,I BD,I CD){ // not sure it is enough
    bool noValidFaces = containsNoTriangle(AB,BC,AC) || containsNoTriangle(AB,BD,AD) || containsNoTriangle(AC,CD,AD)  || containsNoTriangle(BC,CD,BD);
    return noValidFaces;
    /* if (noValidFaces){return true;}
    else {return false;} // forget about the volume
    I AD2 = pow(AD,2);
    I BD2 = pow(BD,2);
    I CD2 = pow(CD,2);
    I X = BD2 + CD2 - pow(BC,2);
    I Y = AD2 + CD2 - pow(AC,2);
    I Z = AD2 + BD2 - pow(AB,2);
    I V = double(4)*AD2*BD2*CD2-AD2*pow(X,2)-BD2*pow(Y,2)-CD2*pow(Z,2)+X*Y*Z;
    bool noVol = (upper(V) <= 0);
    return noVol;    */
}

bool containsNoFM(I rA,I rB,I rC,I rD,I AB,I AC,I AD,I BC,I BD,I CD){
    // True if there are no FM tetrahedra in the set of this t
    return (containsNoTetr(AB,AC,AD,BC,BD,CD) || bad_radius(rA,rB,rC,rD,AB,AC,AD,BC,BD,CD) || lower(supportR(rA,rB,rC,rD,AB,AC,AD,BC,BD,CD)) > upper(r));
}
      
