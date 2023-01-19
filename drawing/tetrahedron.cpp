// to_execute: //g++ tetrahedron.cpp -lglut -lGLU -lGL; ./a.out

#include <GL/glut.h>

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#include "functions_drawing.hpp"

using namespace std;

// Colors
GLfloat BLACK[] = {0, 0, 0,1};
GLfloat BROWN[] = {.2, 0, 0,1};
GLfloat YELLOW[] = { 1.0, 1.0, 0.0, 1.0 };
GLfloat WHITE[] = {1, 1, 1, 1};
GLfloat RED[] = {1, 0, 0, 1};
GLfloat GREEN[] = {0, 1, 0, 1};
GLfloat MAGENTA[] = {1, 0, 1, 1};
GLfloat direction[] = { 1.0, 1.0, 1.0, 0.0 };


void line(float X[3], float Y[3]){
    glBegin(GL_LINES);
    glVertex3f(X[0],X[1],X[2]);
    glVertex3f(Y[0],Y[1],Y[2]);
    glEnd();
}

void tetrahedron(float A[3], float B[3], float C[3], float D[3]){
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, WHITE);
    glColor3f(1.0,1.0,1.0); // blue z
    glLineWidth(2);// line width
    line(A,B);
    line(A,C);
    line(A,D);
    line(B,C);
    line(B,D);
    line(C,D);
    
}

void sphere(float X[3], float rX, GLfloat* col){    
    /*void glutSolidSphere (GLdouble radius , GLint slices , GLint stacks );
      void glutWireSphere(GLdouble radius , GLint slices , GLint stacks );*/
    glLineWidth(.05);// line width
    glPushMatrix();
    
    glLightfv(GL_LIGHT0, GL_DIFFUSE, col);
    
    glTranslatef(X[0], X[1], X[2]);
    glutWireSphere(rX, 20,20);
    glTranslatef(-X[0], -X[1], -X[2]);
    glPopMatrix();
}

void axis(){
    glColor3f(1.0,0.0,0.0); // red x
    glBegin(GL_LINES);
    // x aix   
    glVertex3f(-1.0, 0.0f, 0.0f);
    glVertex3f(1.0, 0.0f, 0.0f);
    glEnd();
    glFlush();
    // y 
    glColor3f(0.0,1.0,0.0); // green y
    glBegin(GL_LINES);
    glVertex3f(0.0, -1.0f, 0.0f);
    glVertex3f(0.0, 1.0f, 0.0f);
    glEnd();
    glFlush();
    
    // z 
    glColor3f(0.0,0.0,1.0); // blue z
    glBegin(GL_LINES);
    glVertex3f(0.0, 0.0f ,-1.0f );
    glVertex3f(0.0, 0.0f ,1.0f );
    glEnd();
    glFlush();
}

float** vertices(float AB, float AC, float AD, float BC, float BD, float CD, float rA, float rB, float rC, float rD){
    if (containsNoFM(rA,rB,rC,rD,AB,AC,AD,BC,BD,CD)){
	cout << "not a valid tetrahedron!!!" << endl;
	return 0;
    }
    float Ax = AD; // d=(0,0,0)
    cout << AD << endl;
    float A[3] = {Ax,0,0}; // a = (ax,0,0)
    float Bx=(pow(AD,2)+pow(BD,2)-pow(AB,2))/(2*AD);
    float By2=pow(BD,2)-pow(Bx,2);
    if (By2<0){
        return 0;
    }
    float By = sqrt(By2);
    float B[3] = {Bx,By,0};
    float Cx=(pow(AD,2)+pow(CD,2)-pow(AC,2))/(2*AD);
    float Cy=-1/2*(pow(AD,2) - pow(AC,2) + pow(BC,2) - pow(Bx,2) - 2*AD*Cx + 2*Bx*Cx - pow(By,2))/By; // wrong...
    float Cz2=pow(CD,2)-pow(Cx,2)-pow(Cy,2);
    if (Cz2<0){
	return 0;
    }
    float Cz=sqrt(Cz2);
    float C[3] = {Cx,Cy,Cz};
    float D[3] = {0,0,0};
    
    float vv[4][3]= {{Ax,0,0},{Bx,By,0},{Cx,Cy,Cz},{0,0,0}};
    float ** v = new float*[4];
    for (int i = 0; i<4; i++){
	v[i] = new float[3];
	for (int j = 0; j<3; j++){
	    v[i][j] = vv[i][j];
	}
    }
    return v;
}


void draw_tetrahedron(float AB, float AC, float AD, float BC, float BD, float CD, float rA, float rB, float rC, float rD, float rot[4], float transl[3]){
    //cout << vertices(1,1,1,1,1,1)[1][0] << endl;
    float ** V = vertices(AB,AC,AD,BC,BD,CD, rA, rB, rC, rD);
    float A[3] = {V[0][0],V[0][1],V[0][2]};
    float B[3] = {V[1][0],V[1][1],V[1][2]};
    float C[3] = {V[2][0],V[2][1],V[2][2]};
    float D[3] = {V[3][0],V[3][1],V[3][2]};
    

    // rotate translate
    glTranslatef(transl[0], transl[1], transl[2]);
    glRotatef(rot[0], rot[1], rot[2], rot[3]); //angle, axis x, y, z   


    glDisable(GL_LIGHTING);        
    glDisable(GL_LIGHT0);          
    glDisable(GL_DEPTH_TEST);   

    tetrahedron(A,B,C,D);
    axis();

    glEnable(GL_LIGHTING);                // so the renderer considers light
    glEnable(GL_LIGHT0);                  // turn LIGHT0 on
    glEnable(GL_DEPTH_TEST);              // so the renderer considers depth*/
    glMaterialfv(GL_FRONT, GL_SPECULAR, WHITE);
    glMaterialf(GL_FRONT, GL_SHININESS, 30);    
    glLightfv(GL_LIGHT0, GL_AMBIENT, BROWN);   
    glLightfv(GL_LIGHT0, GL_SPECULAR, WHITE);
    glLightfv(GL_LIGHT0, GL_POSITION, direction);

    
    sphere(A,rA, YELLOW);
    sphere(B,rB, RED);
    sphere(C,rC, MAGENTA);
    sphere(D,rD, GREEN);

    //  translate rotate
    
    glRotatef(-rot[0], rot[1], rot[2], rot[3]); // go back
    glTranslatef(-transl[0], -transl[1], -transl[2]);


}

void display() {
    float r = 0.414214;
    //--------------------- INPUT -------------------
    // 0.828427,0.828527] [0.828427,0.828527] [0.95052,0.951134] [1.39975,1.40036] [1.2678,1.26842]
    float AB = 2*r, AC = 0.8285,AD = 0.8285,BC = 0.951,BD = 1.3,CD = 1.268;
    //float AB = 2*r, AC = 2*r,AD = 2*r,BC = 2*r,BD = 2*r,CD = 2*r;
    float rA = r, rB = r, rC = r, rD = r;
    float rot1[4] = {300,1,0,0}, rot2[4] = {90,1,0,0}, rot3[4] = {90,1,1,0}, rot4[4] = {300,1,1,0} ;
    float t1[3] = {-1.2,-1.5,0}, t2[3] = {-1.2,1.5,0}, t3[3] = {1.2,1.5,0}, t4[3]{1.2,-1.5,0};
    //cout << supportR(rA,rB,rC,rD,AB,AC,AD,BC,BD,CD) << endl;
    //--------------------- END OF INPUT -------------------
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    
    draw_tetrahedron(AB,AC,AD,BC,BD,CD,rA,rB,rC,rD, rot1, t1);
    
    draw_tetrahedron(AB,AC,AD,BC,BD,CD,rA,rB,rC,rD, rot2, t2);
    
    draw_tetrahedron(AB,AC,AD,BC,BD,CD,rA,rB,rC,rD,rot3, t3);
        
    draw_tetrahedron(AB,AC,AD,BC,BD,CD,rA,rB,rC,rD,rot4, t4);
  
    glPopMatrix();
    glFlush();

}

// We don't want the scene to get distorted when the window size changes, so
// we need a reshape callback.  We'll always maintain a range of -2.5..2.5 in
// the smaller of the width and height for our viewbox, and a range of -10..10
// for the viewbox depth.
void reshape(GLint w, GLint h) {
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  GLfloat aspect = GLfloat(w) / GLfloat(h);
  glLoadIdentity();
  if (w <= h) {
    // width is smaller, so stretch out the height
    glOrtho(-3, 3, -3/aspect, 3/aspect, -10.0, 10.0);
  } else {
    // height is smaller, so stretch out the width
    glOrtho(-3*aspect, 3*aspect, -3, 3, -10.0, 10.0);
  }
}

// Performs application specific initialization.  It defines lighting
// parameters for light source GL_LIGHT0: black for ambient, yellow for
// diffuse, white for specular, and makes it a directional source
// shining along <-1, -1, -1>.  It also sets a couple material properties
// to make cyan colored objects with a fairly low shininess value.  Lighting
// and depth buffer hidden surface removal are enabled here.
void init() {
  GLfloat black[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat yellow[] = { 1.0, 1.0, 0.0, 1.0 };
  GLfloat cyan[] = { 0.0, 1.0, 1.0, 1.0 };
  GLfloat white[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat direction[] = { 1.0, 1.0, 1.0, 0.0 };
    

}

// The usual application statup code.
int main(int argc, char** argv) {
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowPosition(80, 80);
  glutInitWindowSize(800, 600);
  glutCreateWindow("Tetrahedron");
  glutReshapeFunc(reshape);
  glutDisplayFunc(display);
  init();
  glutMainLoop();
}
