CC = gcc
CCFLAGS = -Wall -Werror -pedantic -std=c99 -g
METIS_INC = -I/homes/kt208/metis/metis_install/include/
METIS_LDINC = -L/homes/kt208/metis/metis_install/lib/
LDFLAGS = -lm -lmetis

mesh: mesh_part.c airfoil_utils.h 
	$(CC) -o $@ $(CCFLAGS) $(METIS_INC) $(METIS_LDINC) $< $(LDFLAGS)

airfoil: airfoil.cpp
	g++ -o $@ -O3 -lm $<

%.svg: %.dot
	./dot2svg $^

%.dot: mesh
	./$^

graphs: meshColoured.svg

clean:
	rm -f mesh
	rm -f airfoil
	rm -f *.svg
