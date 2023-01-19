GLFLAGS =  -lglut -lGLU -lGL

all: draw

draw: tetrahedron.o functions_drawing.o
	g++ -o draw   tetrahedron.o functions_drawing.o  ${GLFLAGS} 

tetrahedron.o: tetrahedron.cpp 
	g++ ${GLFLAGS} -c tetrahedron.cpp 

functions_drawing.o: functions_drawing.cpp
	g++ -c functions_drawing.cpp

clean:
	rm -f draw tetrahedron.o functions_drawing.o
