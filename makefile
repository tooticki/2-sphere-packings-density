all: program

program: main.o functions.o
	g++ -o test main.o functions.o

main.o: main.cpp
	g++ -c -ggdb3 main.cpp

functions.o: functions.cpp
	g++ -c -ggdb3 functions.cpp

clean:
	rm -f test main.o functions.o
