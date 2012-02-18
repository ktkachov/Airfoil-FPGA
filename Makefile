all: airfoil.cpp
	g++ -O3 -o airfoil airfoil.cpp

mesh: mesh_part.cpp
	gcc -o mesh -I/homes/kt208/metis/metis_install/include -L/homes/kt208/metis/metis_install/lib/ mesh_part.cpp -lm -lmetis

clean:
	rm -f mesh
	rm -f airfoil
