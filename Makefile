CC = gcc
CCFLAGS = -Wall -Werror -pedantic -std=c99 -g
CCINC = -I/homes/kt208/metis/metis_install/include/
LDINC = -L/homes/kt208/metis/metis_install/lib/
LDFLAGS = -lm -lmetis

all: mesh_part.c hash_map.h 
	$(CC) -o mesh $(CCFLAGS) $(CCINC) $(LDINC) $< $(LDFLAGS)

airfoil: airfoil.cpp
	g++ -o airfoil -O3 -lm $<

clean:
	rm -f mesh
	rm -f airfoil
