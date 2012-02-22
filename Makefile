CC = gcc
CCFLAGS = -g -Wall -Werror -pedantic -std=c99
CCINC = -I/homes/kt208/metis/metis_install/include/
LDINC = -L/homes/kt208/metis/metis_install/lib/
LDFLAGS = -lm -lmetis

airfoil: airfoil.cpp
	g++ -O3 -o airfoil airfoil.cpp

mesh: mesh_part.c
	$(CC) -o mesh $(CCFLAGS) $(CCINC) $(LDINC) mesh_part.c $(LDFLAGS)

clean:
	rm -f mesh
	rm -f airfoil
