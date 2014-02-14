#!/bin/make -f

#LDFLAGS=
#LDFLAGS=-L../src/.libs/ -static
LDFLAGS=-L../src/.libs/ 
#CFLAGS=-c -Wall -Wno-deprecated -I../include/
#CFLAGS=-c -O3 -fno-inline -ffast-math -Wall -Wno-deprecated -I../include/
#CFLAGS=-c -g -O3 -ffast-math -Wall -Wno-deprecated -I../include/
#CFLAGS= -g -pg -O3 -flto --fast-math -fno-inline  -Wall -Wno-deprecated -I../include/
CFLAGS=-O3 -g -fopenmp -march=native -msse3 -flto --fast-math -Wall -Wno-deprecated -I../include/
#CFLAGS=-g -Wall -Wno-deprecated -I../include/
#CFLAGGS=-O3 -g -ipo -msse4.1 -march=native  -Wno-deprecated -Wall -I../include/

INCLUDES=-I../include/ 
#CC=LC_ALL=C icc 
#CC=icpc
CC= g++
LIBS=-lm -lstdc++ -lgsl -lgslcblas 


all: McPoly

McPoly:  polyhedra.o cube.o octahedron.o cuboid.o truncated_cube.o particles.o move.o Translate_Brownian.o Rotate_2.o box.o fileio.o cell_list.o cell.o main3.o order_parameter_2.o cluster_1.o Vol_move.o alg_struct.o Collision_Test.o tis_window.o
	$(CC) $(CFLAGS) polyhedra.o cube.o octahedron.o particles.o cuboid.o truncated_cube.o move.o Translate_Brownian.o Rotate_2.o box.o fileio.o cell_list.o cell.o main3.o order_parameter_2.o cluster_1.o Vol_move.o alg_struct.o Collision_Test.o tis_window.o $(INCLUDES) $(LDFLAGS) $(LIBS) -o McPoly 

cube.o: cube.cpp
	$(CC) $(CFLAGS) -c cube.cpp  

polyhedra.o: polyhedra.cpp
	$(CC) $(CFLAGS) -c polyhedra.cpp

particles.o: particles.cpp
	$(CC) $(CFLAGS) -c particles.cpp  

cuboid.o: cuboid.cpp
	$(CC) $(CFLAGS) -c cuboid.cpp 

move.o: move.cpp
	$(CC) $(CFLAGS) -c move.cpp

Rotate_2.o: Rotate_2.cpp
	$(CC) $(CFLAGS) -c Rotate_2.cpp 
 

Translate_Brownian.o: Translate_Brownian.cpp
	$(CC) $(CFLAGS) -c Translate_Brownian.cpp 

box.o: box.cpp
	$(CC) $(CFLAGS) -c box.cpp  	

fileio.o: fileio.cpp
	$(CC) $(CFLAGS) -c fileio.cpp 

cell_list.o: cell_list.cpp
	$(CC) $(CFLAGS) -c cell_list.cpp 

cell.o: cell.cpp
	$(CC) $(CFLAGS) -c cell.cpp 

main3.o: main3.cpp
	$(CC) $(CFLAGS) -c main3.cpp 

order_parameter_2.o: order_parameter_2.cpp
	$(CC) $(CFLAGS) -c order_parameter_2.cpp   	

cluster_1.o: cluster_1.cpp
	$(CC) $(CFLAGS) -c cluster_1.cpp

Vol_move.o: Vol_move.cpp
	$(CC) $(CFLAGS) -c Vol_move.cpp   

alg_struct.o: alg_struct.cpp
	$(CC) $(CFLAGS) -c alg_struct.cpp   

Collision_Test.o: Collision_Test.cpp
	$(CC) $(CFLAGS) -c Collision_Test.cpp   

tis_window.o: tis_window.cpp
	$(CC) $(CFLAGS) -c tis_window.cpp

octahedron.o: octahedron.cpp
	$(CC) $(CFLAGS) -c octahedron.cpp

truncated_cube.o: truncated_cube.cpp
	$(CC) $(CFLAGS) -c truncated_cube.cpp


clean:
	rm -rf *o McPoly


