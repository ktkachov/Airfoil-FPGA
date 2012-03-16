CC = gcc
CCFLAGS = -Wall -Werror -pedantic -std=c99 -g
CCINC = -I/homes/kt208/metis/metis_install/include/
LDINC = -L/homes/kt208/metis/metis_install/lib/
LDFLAGS = -lm -lmetis

mesh: mesh_part.c airfoil_utils.h 
	$(CC) -o $@ $(CCFLAGS) $(CCINC) $(LDINC) $< $(LDFLAGS)

airfoil: airfoil.cpp
	g++ -o $@ -O3 -lm $<

svg:
	./dot2svg *.dot

clean:
	rm -f mesh
	rm -f airfoil
	rm -f *.svg
