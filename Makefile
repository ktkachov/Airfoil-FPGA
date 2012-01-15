all: airfoil.cpp
	g++ -O3 -o airfoil airfoil.cpp

clean:
	rm airfoil
