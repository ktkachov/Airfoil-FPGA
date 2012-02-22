airfoil: airfoil.cpp
	g++ -O3 -o airfoil airfoil.cpp

mesh: mesh_part.c
	gcc -o mesh -g -std=c99 -I/homes/kt208/metis/metis_install/include -L/homes/kt208/metis/metis_install/lib/ mesh_part.c -lm -lmetis

clean:
	rm -f mesh
	rm -f airfoil
