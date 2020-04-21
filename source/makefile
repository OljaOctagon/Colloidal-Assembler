#!/bin/make -f
CC=g++
LIBS=-lgsl -lgslcblas 
CFLAGS=-O3 -g

all: McPoly

McPoly:  Three_Particle_Move.o Cluster_Move.o Grand_Canonical.o polyhedra.o cube.o octahedron.o cuboid.o truncated_cube.o rhombohedron.o rectangle.o hexbipyramid.o particles.o pmove.o translate.o rotate_3D.o box.o fileio.o cell_list.o cell.o main.o order_parameter.o cluster.o Vol_move.o alg_struct.o Collision_Test.o tis_window.o parallelepiped.o rotate_2D.o gear.o
	$(CC) $(CFLAGS) Three_Particle_Move.o Cluster_Move.o Grand_Canonical.o polyhedra.o cube.o octahedron.o particles.o cuboid.o truncated_cube.o rhombohedron.o rectangle.o hexbipyramid.o pmove.o translate.o rotate_3D.o box.o fileio.o cell_list.o cell.o main.o order_parameter.o cluster.o Vol_move.o alg_struct.o Collision_Test.o tis_window.o parallelepiped.o rotate_2D.o gear.o $(LIBS) -o McPoly 

gear.o: gear.cpp
	$(CC) $(CFLAGS) -c gear.cpp

cube.o: cube.cpp
	$(CC) $(CFLAGS) -c cube.cpp  

polyhedra.o: polyhedra.cpp
	$(CC) $(CFLAGS) -c polyhedra.cpp

particles.o: particles.cpp
	$(CC) $(CFLAGS) -c particles.cpp  

cuboid.o: cuboid.cpp
	$(CC) $(CFLAGS) -c cuboid.cpp 

pmove.o: pmove.cpp
	$(CC) $(CFLAGS) -c pmove.cpp

rotate_3D.o: rotate_3D.cpp
	$(CC) $(CFLAGS) -c rotate_3D.cpp 

rotate_2D.o: rotate_2D.cpp
	$(CC) $(CFLAGS) -c rotate_2D.cpp

translate.o: translate.cpp
	$(CC) $(CFLAGS) -c translate.cpp 

box.o: box.cpp
	$(CC) $(CFLAGS) -c box.cpp

fileio.o: fileio.cpp
	$(CC) $(CFLAGS) -c fileio.cpp 

cell_list.o: cell_list.cpp
	$(CC) $(CFLAGS) -c cell_list.cpp 

cell.o: cell.cpp
	$(CC) $(CFLAGS) -c cell.cpp 

main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp 

order_parameter.o: order_parameter.cpp
	$(CC) $(CFLAGS) -c order_parameter.cpp   	

cluster.o: cluster.cpp
	$(CC) $(CFLAGS) -c cluster.cpp

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

rhombohedron.o: rhombohedron.cpp
	$(CC) $(CFLAGS) -c rhombohedron.cpp

rectangle.o: rectangle.cpp
	$(CC) $(CFLAGS) -c rectangle.cpp

hexbipyramid.o: hexbipyramid.cpp
	$(CC) $(CFLAGS) -c hexbipyramid.cpp

parallelepiped.o: parallelepiped.cpp
	$(CC) $(CFLAGS) -c parallelepiped.cpp

Grand_Canonical.o: Grand_Canonical.cpp
	$(CC) $(CFLAGS) -c Grand_Canonical.cpp 

Cluster_Move.o: Cluster_Move.cpp
	$(CC) $(CFLAGS) -c Cluster_Move.cpp

Three_Particle_Move.o: Three_Particle_Move.cpp
	$(CC) $(CFLAGS) -c Three_Particle_Move.cpp

clean:
	rm -rf *o McPoly
